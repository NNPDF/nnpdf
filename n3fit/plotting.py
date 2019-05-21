#!/usr/bin/env python3

"""
    Script to plot the hyperopt scans
"""

import os
import json
import warnings
import re

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

warnings.simplefilter(action="ignore", category=FutureWarning)
regex_op = re.compile(r"[^\w^\.]+")
regex_not_op = re.compile(r"[\w\.]+")

keywords = {
        'id' : 'iteration',

        'optimizer' : 'optimizer',
        'lr' : 'learning_rate',
        'initializer' : 'initializer',
        'dropout' : 'dropout',

        'nodes' : 'nodes_per_layer',
        'max_layers' : 4,
        'nl' : 'number_of_layers',
        'activation' : 'activation_per_layer',
        'architecture' : 'layer_type',

        'epochs' : 'epochs',
        'stp' : 'stopping_patience',
        'ste' : 'stopping_epochs',
        'p_ini' : 'pos_initial',
        'p_mul' : 'pos_multiplier',

        'good' : 'good',
        'vl' : 'validation_loss',
        'tl' : 'testing_loss',
        'loss' : 'loss',
        }

# 0 = normal scatter plot, 1 = violin, 2 = log
plotting_keys = [
        (keywords['id'], 0),
        (keywords['optimizer'], 1),
        (keywords['lr'], 2, [2e-4, 4e-1]),
        (keywords['initializer'], 1),

        (keywords['epochs'], 0),
        (keywords['ste'], 0),
        (keywords['stp'], 0),
        (keywords['p_mul'], 0),

        (keywords['nl'], 1),
        (keywords['activation'], 1),
        (keywords['dropout'], 0),
        (keywords['tl'], 0),

        ]

# for i in range(4):
#     plotting_keys.append( ('layer_{0}'.format(i+1),0) )


# Family parsing
def optimizer_parse(trial, dict_out):
    """
    This function parses the family of parameters which regards the optimizer:

    optimizer
    learning_rate (when it exists)
    initializer
    dropout
    """
    opt = trial["misc"]["space_vals"][keywords["optimizer"]]
    if isinstance(opt, dict):
        # if it is a dict the optimizer includes extra information
        name = opt[keywords["optimizer"]]
        lr = opt[keywords["lr"]]
    else:
        name = opt
        lr = None
    ini = trial["misc"]["space_vals"][keywords["initializer"]]
    dropout_rate = trial["misc"]["space_vals"][keywords["dropout"]]

    dict_out[keywords["optimizer"]] = name
    dict_out[keywords["lr"]] = lr
    dict_out[keywords["initializer"]] = ini
    dict_out[keywords["dropout"]] = dropout_rate


def architecture_parser(trial, dict_out):
    """
    This function parses the family of parameters which regards the architecture of the NN

    number_of_layers
    activation_per_layer
    nodes_per_layer
    l1, l2, l3, l4... max_layers
    layer_type
    """
    nodes_per_layer = trial["misc"]["space_vals"][keywords["nodes"]]
    nl = len(nodes_per_layer) - 1
    activation_name = trial["misc"]["space_vals"][keywords["activation"]][0]
    architecture = trial["misc"]["space_vals"][keywords["architecture"]]

    dict_out[keywords["nodes"]] = nodes_per_layer
    dict_out[keywords["nl"]] = nl
    dict_out[keywords["activation"]] = activation_name
    dict_out[keywords["architecture"]] = architecture

    # In principle we might be checking any number of layers,
    # so it is important to make sure that the maximum number of layers
    # (which will be used at plotting time) is correct
    if nl > keywords["max_layers"]:
        keywords["max_layers"] = nl

    for i in range(keywords["max_layers"]):
        dict_out["layer_{0}".format(i + 1)] = None

    for i, nodes in enumerate(nodes_per_layer[:-1]):
        dict_out["layer_{0}".format(i + 1)] = nodes


def stopping_parser(trial, dict_out):
    """
    Parses the family of parameters that define the stopping

    epochs
    stopping_patience
    pos_initial
    pos_multiplier
    """

    epochs = trial["misc"]["space_vals"][keywords["epochs"]]
    patience = trial["misc"]["space_vals"][keywords["stp"]]
    stop_ep = patience * epochs
    positivity_initial = trial["misc"]["space_vals"][keywords["p_ini"]]
    positivity_multiplier = trial["misc"]["space_vals"][keywords["p_mul"]]

    dict_out[keywords["epochs"]] = epochs
    dict_out[keywords["stp"]] = patience
    dict_out[keywords["ste"]] = stop_ep
    dict_out[keywords["p_ini"]] = positivity_initial
    dict_out[keywords["p_mul"]] = positivity_multiplier


def stat_parser(trial, dict_out, threshold=1e3, val_f=0.5, **kwargs):
    """
    Parser the stats information of the trial, the losses, the failures, etc
    """
    test_f = 1.0 - val_f
    validation_loss = trial["result"][keywords["vl"]]
    testing_loss = trial["result"][keywords["tl"]]
    loss = validation_loss * val_f + testing_loss * test_f

    if testing_loss < threshold and validation_loss < threshold:
        good_result = True
    else:
        good_result = False

    # Was this a ok run?
    if trial["result"]["status"] != "ok":
        ok = False
    else:
        ok = good_result

    if not good_result:
        # Set bad results to be an order of magnitude above everything else
        loss = threshold * 10

    dict_out[keywords["good"]] = ok
    dict_out[keywords["loss"]] = loss
    dict_out[keywords["vl"]] = validation_loss
    dict_out[keywords["tl"]] = testing_loss


def filter_by_string(data_dict, filter_string):
    """
    Receives a data_dict (a parsed trial) and a filter string, returns True if the trial passes the filter

    filter string must have the format: key[]string
    where [] can be any of !=, =, >, <
    """
    match = regex_not_op.findall(filter_string)
    if len(match) < 2:
        raise Exception("Filter str is not correct: {0}".format(filter_string))
    filter_key = match[0]
    filter_val = match[1]

    if filter_key in keywords.keys():
        filter_key = keywords[filter_key]

    val_check = data_dict[filter_key]
    if val_check is None:  # NaN means it does not apply
        return True

    operator = regex_op.findall(filter_string)[0]

    if operator == "=":
        operator = "=="
    operators = ["!=", "==", ">", "<"]
    if operator not in operators:
        raise Exception("Filter string not valid, operator not recognized {0}".format(filter_string))

    # This I know it is not ok:
    if isinstance(val_check, str) and isinstance(filter_val, str):
        check_str = '"{0}"{1}"{2}"'
    else:
        check_str = "{0}{1}{2}"
    try:
        return eval(check_str.format(val_check, operator, filter_val))
    except:
        return False


# General parsing
def trial_parser(trial, filter_me=None, **kwargs):
    """
    Trials are very convoluted object, very branched inside
    The goal of this function is to separate said branching so we can create hierarchies
    """
    # Is this a true trial
    if trial["state"] != 2:
        return None

    data_dict = {}
    # Parse all information in easy to use dictionaries
    optimizer_parse(trial, data_dict)
    architecture_parser(trial, data_dict)
    stopping_parser(trial, data_dict)
    stat_parser(trial, data_dict, **kwargs)

    # If a filter is set, check whether the trial pass the filter
    if filter_me:
        for filter_str in filter_me:
            if not filter_by_string(data_dict, filter_str):
                return None

    return data_dict


def filter_trials(input_trials, **kwargs):
    """
    Applies the search algorithm to all the different trials to find the true best model
    """
    # First step is to fill in two different
    failed_trials = []
    ok_trials = []

    full_info = []
    for tid, trial in enumerate(input_trials):
        trial_dict = trial_parser(trial, **kwargs)
        if trial_dict is None:
            continue
        trial_dict[keywords["id"]] = tid
        full_info.append(trial_dict)

    return full_info


def order_axis(df, bestdf, key):
    """
    Helper function for ordering the axis and make sure the best is always first
    """
    best_x_lst = bestdf.get(key).tolist()
    ordering = set(df.get(key).tolist())
    ordering.remove(best_x_lst[0])
    ordering_true = best_x_lst + list(ordering)
    best_x = np.array([str(best_x_lst[0])])
    return ordering_true, best_x


def plot_scans(df, best_df, outfile):
    """
    This function plots all trials in a nice multiplot
    """
    # Some printouts
    print("All trials:")
    print(df)

    print("Best setup:")
    with pd.option_context("display.max_rows", None, "display.max_columns", None):
        print(best_df)

    # Append extra keys for the number of possible layers
    for i in range(keywords["max_layers"]):
        plotting_keys.append(("layer_{0}".format(i + 1), 4))

    # Create a grid of plots
    nplots = len(plotting_keys)
    _, axs = plt.subplots(4, nplots // 4, sharey=True, figsize=(30, 30))

    # Set the quantity we will be plotting in the y axis
    loss_k = keywords["loss"]

    for ax, key_tuple in zip(axs.reshape(-1), plotting_keys):
        key = key_tuple[0]
        mode = key_tuple[1]

        if mode == 0 or mode == 2:  # normal scatter plot
            ax = sns.scatterplot(x=key, y=loss_k, data=df, ax=ax)
            best_x = best_df.get(key)
            if mode == 2:
                ax.set_xscale("log")
                ax.set_xlim(key_tuple[2])
        elif mode == 1:
            ordering_true, best_x = order_axis(df, best_df, key=key)
            ax = sns.violinplot(x=key, y=loss_k, data=df, ax=ax, palette="Set2", cut=0.0, order=ordering_true)
            ax = sns.stripplot(x=key, y=loss_k, data=df, ax=ax, color="gray", order=ordering_true, alpha=0.4)
        elif mode == 4:
            # The nodes per layer is a special case because we want to know to
            # which total number of layers they correspond
            ordering_true, _ = order_axis(df, best_df, key=keywords["nl"])
            best_x = best_df.get(key)
            for l in ordering_true:
                plot_data = df[df[keywords["nl"]] == l]
                label = "layers = {0}".format(l)
                ax = sns.scatterplot(x=key, y=loss_k, data=plot_data, ax=ax, label=label)
            ax.legend()

        # Finally plot the "best" one, which will be first
        ax = sns.scatterplot(x=best_x, y=best_df.get(loss_k), ax=ax, color="orange", marker="s")
        ax.set_ylabel("Loss")
        ax.set_xlabel(key)

    plt.savefig(f"{outfile}", bbox_inches="tight")


def draw_histogram(dataframe, key, plotable="loss", columns=2):
    values = set(dataframe[key])
    nv = len(values)
    rows = np.ceil(nv / columns)

    for i, value in enumerate(values):
        plt.subplot(rows, columns, i + 1)
        data = dataframe[dataframe[key] == value]
        plt.hist(data[plotable], 20)
        plt.xlabel(value)
    plt.show()


def print_stats(
    dataframe,
    stats_keys=["optimizer", "nl", "initializer", "activation"],
    n_combinations=2,
    verbose=True,
    return_early=False,
    how_many_to_kill=2,
    remove_by_key=False,
    save_n_best=10,
    **kwargs,
):
    """
    Print statistics on failure rates, averages and so
    """

    # Ok, at this point I have a dataframe with all the content from the json.
    # if filters or changes on the threshold for the loss are active, the dataframe would
    # already be filtered and all information ready

    def read_key(key):
        if key in keywords.keys():
            key_name = keywords[key]
        else:
            key_name = key
        return key_name

    def kclean(l):
        return [read_key(i) for i in l]

    def remove_key(dictionary, key):
        new_d = dictionary.copy()
        new_d.pop(key)
        return new_d

    def gen_key_info(full_dataframe, keys):
        key_info = {}
        for key_name in keys:
            all_possible_values = full_dataframe[key_name]
            remove_nans = filter(lambda x: x == x, all_possible_values)
            remove_duplicates = sorted(list(set(remove_nans)))
            if remove_duplicates:
                key_info[key_name] = remove_duplicates
        return key_info

    def get_slice(full_dataframe, key, value):
        return full_dataframe[full_dataframe[key] == value]

    from HyperAlgorithm import HyperAlgorithm

    keys_lv_0 = kclean(["optimizer", "number_of_layers"])
    keys_lv_1 = kclean(["epochs", "stopping_patience", "activation"])
    keys_lv_2 = kclean(["p_mul", "p_ini", "learning_rate"])

    keys_for_removal = keys_lv_0 + keys_lv_1

    key_info_full = gen_key_info(dataframe, keys_for_removal)
    hypalg = HyperAlgorithm(key_info=key_info_full)
    trim_dataframe = hypalg.remove_failing(dataframe)
    key_info_trim = gen_key_info(trim_dataframe, keys_for_removal)
    hypalg = HyperAlgorithm(key_info=key_info_trim)

    decision_tree = {}

    key_0 = keys_lv_0[0]
    options_0 = key_info_trim[key_0]
    for val in options_0:
        print("For {0}".format(val))
        full_slice = get_slice(trim_dataframe, key_0, val)
        trim_slice = hypalg.remove_failing(full_slice)
        for key in keys_for_removal[1:]:
            trim_slice = hypalg.remove_failing(trim_slice, key=key, fail_threshold=0.0)
        # At this point I've already removed the worst worst ones, now we can create the "decision_tree"
        tree_key_info = gen_key_info(trim_slice, keys_lv_0[0:])
        decision_tree[val] = (tree_key_info, trim_slice)

    print("\n\n")
    verbose = False

    final_decision_tree = {}
    key_1 = keys_lv_0[1]
    for val_0, info in decision_tree.items():
        keys = info[0]
        f_slice = info[1]
        options_1 = keys[key_1]
        if verbose:
            print("For {0}, accepted options are: {1}={2}".format(val_0, key_1, options_1))
        final_decision_tree[val_0] = {}
        for val in options_1:
            print(val)
            trim_slice = get_slice(f_slice, key_1, val)
            for key_name in keys_for_removal[2:]:
                trim_slice = hypalg.remove_failing(trim_slice, key=key_name, fail_threshold=0.0, verbose=False)
            # Save the final slice for this branch of the tree
            final_decision_tree[val_0][val] = trim_slice

    # Now for each branch let's look at the rewards
    final_answer = {}
    for val_0, item_0 in final_decision_tree.items():
        if verbose:
            print("For {0}".format(val_0))
        final_answer[val_0] = {}
        for val_1, item_1 in item_0.items():
            if verbose:
                print("\n\nFor {0}".format(val_1))
            # Now remove half of the combinations but only looking at the level_1 keys
            key_info_lv1 = gen_key_info(item_1, keys_lv_1)
            falg = HyperAlgorithm(key_info=key_info_lv1)
            n1 = len(key_info_lv1)
            trim_slice = falg.remove_failing(item_1, n_combinations=n1, fail_threshold="median", verbose=False)
            # Now do the same with the level 2, remove the worst half
            key_info_lv2 = gen_key_info(trim_slice, keys_lv_2)
            falg = HyperAlgorithm(key_info=key_info_lv2)
            n2 = len(key_info_lv2)
            final_slice = falg.remove_failing(trim_slice, n_combinations=n2, fail_threshold="median", verbose=False)
            # Now print the ranges
            print_ranges = keys_lv_1 + keys_lv_2
            print(
                """
For {0} = {1}
    {2} = {3}""".format(
                    key_0, val_0, key_1, val_1
                )
            )
            palg = HyperAlgorithm(key_info=gen_key_info(final_slice, print_ranges))
            palg.print_ranges(final_slice, dataframe_original=item_1)

    return

    # key_priority_raw = ["architecture", "optimizer", "number_of_layers", "epochs", "activation"]
    key_priority_raw = [
        "optimizer",
        "epochs",
        "number_of_layers",
        "learning_rate",
        "activation",
        "layer_1",
        "layer_2",
        "layer_3",
    ]
    key_priority = [read_key(i) for i in key_priority_raw]
    # Now, for each key we need to get all the possible options

    key_info = gen_key_info(dataframe)
    # Now that we have all the options we need to go in priority order

    # 1. Remove the failing entries (i.e., reward == failing reward)
    hypalg = HyperAlgorithm(key_info=key_info)
    trim_dataframe = hypalg.remove_failing(dataframe)
    # 1.b Regenerate key_info
    key_info = gen_key_info(trim_dataframe)

    # 2. Now fix the optimizer and continue
    key_0 = key_priority[0]  # optimizer
    options_0 = key_info[key_0]
    key_info_0 = remove_key(key_info, key_0)
    hypalg = HyperAlgorithm(key_info=key_info_0)
    for val in options_0:
        print("\n\n")
        print(val)
        f_slice = trim_dataframe[trim_dataframe[key_0] == val]
        # First remove all failing keys
        trim_slice = hypalg.remove_failing(f_slice)
        for key_name in key_priority[1:]:
            print(">")
            print(key_name)
            trim_slice = hypalg.remove_failing(trim_slice, key=key_name, fail_threshold=0.0)
        hypalg.print_ranges(trim_slice, dataframe_original=f_slice)

    set_trace()

    # Create a key_info dict with one entry removed for each key
    key_infos = []
    last_k = key_info
    for key_name in key_priority:
        last_k = remove_key(last_k, key_name)
        key_infos.append(last_k)

    key_0 = key_priority[0]
    options_0 = key_info[key_0]
    slices_0 = []
    key_info_0 = {key_priority[1]: key_info[key_priority[1]]}
    hyperalg_0 = HyperAlgorithm(key_info=key_info)
    hyperalg_0.check_wrapper_weaklings_killer(dataframe, 1)
    for option_0 in options_0:
        slices_0.append(dataframe[dataframe[key_0] == option_0])
    for slice_0 in slices_0:
        # check the rewards, if some are particularly bad, kill them
        _ = hyperalg_0.check_wrapper_weaklings_killer(slice_0, 1)
    # entonce smatamos primero los epochs por ejemplo
    # con los epochs malos eliminados, miramos el numbeor de layers
    # y eliminamos los dos peores (nos quedamos con la mitad)
    # luego para cada uno nos quedamos con el set de numero de nodos que tiene mejor pinta
    # y asi
    set_trace()

    def print_me(string):
        if verbose:
            print(string)

    def kill_val(key, name):
        dataframe.loc[dataframe[key] == name, keywords["good"]] = False

    def process_slice(df_slice):
        """ Function to process a slice into a dictionary with interesting stats """
        md = {}
        good = df_slice[df_slice[keywords["good"]]]
        md["n_total"] = len(df_slice)
        md["n_good"] = len(good)
        md["n_failed"] = md["n_total"] - md["n_good"]
        if md["n_total"] == 0:
            f_rate = None
        else:
            f_rate = md["n_failed"] / md["n_total"] * 100
        md["f_rate"] = f_rate

        # Compute for each the best loss (among the good ones) and the std dev
        good_losses = good[keywords["loss"]]
        md["best_loss"] = good_losses.min()
        md["std"] = good_losses.std()
        return md

    ###############################

    # Parse the keys into a dict of { 'keyname' : [possible, values] }
    key_info = {}
    for key in stats_keys:
        if key in keywords.keys():
            key_name = keywords[key]
        else:
            key_name = key
        possible_values = list(
                set(
                    filter( lambda x: x == x, dataframe[key_name] )
                )
            )
        key_info[key_name] = possible_values

    from HyperAlgorithm import HyperAlgorithm

    hyperalg = HyperAlgorithm(
        key_info=key_info, fail_threshold=65, how_many_to_kill=how_many_to_kill, remove_by_key=remove_by_key
    )

    # First run over all keys and print_me the information for them to see the general picture

    # Perform a genocide for every reward below the threshold for 1 combination and then 2 the same with combinations of 2
    new_dataframe = hyperalg.check_wrapper_weaklings_killer(dataframe, 1)
    new_dataframe = hyperalg.check_wrapper_weaklings_killer(new_dataframe, 2)
    set_trace()
    new_dataframe = hyperalg.check_wrapper_weaklings_killer(new_dataframe, 1)
    new_dataframe = hyperalg.check_wrapper_weaklings_killer(new_dataframe, 2)
    set_trace()
    new_dataframe = hyperalg.check_wrapper_weaklings_killer(new_dataframe, 3)

    #     for i in range(1, n_combinations):
    #         new_dataframe = hyperalg.check_wrapper(new_dataframe, i, kill = True)
    # For the last one, only the n_best surviv
    new_dataframe = hyperalg.check_wrapper(new_dataframe, n_combinations, kill=True, save_n_best=save_n_best)

    # Finally print the survivors
    hyperalg.print_ranges(new_dataframe, dataframe_original=dataframe)


# Wrapper
def generate_scan_report_from_file(
    replica_path, json_name="tries.json", include_failures=False, histogram=None, stats=False, **kwargs
):
    """
    Read json file and generate a the corresponding scan report
    """
    filename = "{0}/{1}".format(replica_path, json_name)

    # Open the file and reads it as a json
    with open(filename, "r") as jlist:
        input_trials = json.load(jlist)

    # Read all trials and create a list of dictionaries
    trials_dictionary = filter_trials(input_trials, **kwargs)

    # Make the list of dictionaries into a dataframe
    dataframe_raw = pd.DataFrame(trials_dictionary)

    # If we just want to print stats, do that and exit
    if stats:
        print_stats(dataframe_raw, **kwargs)
        return

    # Now select whether we want only good trials or rather all of them
    if include_failures:
        dataframe = dataframe_raw
    else:
        dataframe = dataframe_raw[dataframe_raw["good"]]

    # Now select the best one
    best_idx = dataframe.loss.idxmin()
    best_trial = dataframe.ix[best_idx]
    # Select it in  this a bit more complicatted way for the plotting code to work ????
    best_trial = dataframe[dataframe[keywords["id"]] == best_trial.get(keywords["id"])]

    if histogram:
        draw_histogram(dataframe, histogram)

    # Plot the plot
    fileout = f"{replica_path}/scan.pdf"
    plot_scans(dataframe, best_trial, fileout)
    print("Plot saved at {0}".format(fileout))


if __name__ == "__main__":
    from argparse import ArgumentParser
    import glob

    parser = ArgumentParser()
    parser.add_argument("folder", help = "Fit folder")
    parser.add_argument("-v", "--val_multiplier", help = 'Fraction to weight the validation loss with (test_multiplier = 1-val_multiplier)',
            type = float, default = 0.5)
    parser.add_argument("-if", "--include_failures", help = 'Flag to include failed runs in  the plots', action = 'store_true')
    parser.add_argument("-t", "--threshold", help = 'Value of the loss function from which to consider a run to have failed',
            type = float, default = 1e3)
    parser.add_argument("-dh", "--histogram", help = 'Write a histogram for a given key')
    parser.add_argument("-f", "--filter", help = "Add the filter key=value to the dataframe", nargs = '+')
    parser.add_argument("-s", "--stats", help = "Print stats on the run instead of plotting", action = 'store_true')
    parser.add_argument("-k", "--stats_keys", help = "When using toghether with -s, chooses which keys to look at",
            default = ['optimizer', 'number_of_layers', 'initializer', 'layer_1', 'layer_2', 'layer_3', 'layer_4', 'learning_rate', 'epochs'],
            nargs = '+')
    parser.add_argument("--ncombinations", help = "When using together with -s, chooses how many keys to combine", type = int, default = 2)
    parser.add_argument("--n_to_remove", help = "Number of worst values (by reward) to remove", type = int, default = 2)
    parser.add_argument("--remove_by_key", help = "Remove by key instead of by the worst of all combinations", action = "store_true")
    parser.add_argument("--save_n_best", help = "Prints the n best models as a range", type = int, default = 10)
    args = parser.parse_args()

    # Search for all json for all replicas
    search_str = "{0}/nnfit/replica_*/tries.json".format(args.folder)
    all_json = glob.glob(search_str)

    for json_path in all_json:
        replica_path = os.path.dirname(json_path)
        generate_scan_report_from_file(
            replica_path,
            val_f=args.val_multiplier,
            include_failures=args.include_failures,
            threshold=args.threshold,
            histogram=args.histogram,
            filter_me=args.filter,
            stats=args.stats,
            stats_keys=args.stats_keys,
            n_combinations=args.ncombinations,
            how_many_to_kill=args.n_to_remove,
            remove_by_key=args.remove_by_key,
            save_n_best=args.save_n_best,
        )
