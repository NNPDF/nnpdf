#!/usr/bin/env python3

"""
    This module contains the necessary functions for generating
    the hyperopt plots.

    In the most basic form it simply generates a plot with all the trials
    with the lowest loss configuration drawn in a different color.

    Receives as input a fit folder containing at least one replica.
    If there are more than one replicas in the folder it will run iteratively
    through all the json unless the option --combine is given

    It also includes several options to massage the results:

    # Options:
    -f, --filter: expresions of the form key [operator] value
              the allowed operators are: =, <, >, !=
    -v: change the value of the validation multiplier (by default 0.5)
        the test multiplier is just (1.0 - v)
    -t: if the loss is above this value, the configuration is considered as failure
    -if: if given, the flags include failing runs
    -c: if given, combine all the replicas as if it were one big json file
"""

from argparse import ArgumentParser
import glob
import json
import os
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

from n3fit.hyper_optimization.hyper_algorithm import autofilter_dataframe

regex_op = re.compile(r"[^\w^\.]+")
regex_not_op = re.compile(r"[\w\.]+")

# ...
keywords = {
    "id": "iteration",
    "optimizer": "optimizer",
    "lr": "learning_rate",
    "initializer": "initializer",
    "dropout": "dropout",
    "nodes": "nodes_per_layer",
    "max_layers": 4,
    "nl": "number_of_layers",
    "activation": "activation_per_layer",
    "architecture": "layer_type",
    "epochs": "epochs",
    "stp": "stopping_patience",
    "ste": "stopping_epochs",
    "p_ini": "pos_initial",
    "p_mul": "pos_multiplier",
    "good": "good",
    "vl": "validation_loss",
    "tl": "testing_loss",
    "loss": "loss",
}

# 0 = normal scatter plot, 1 = violin, 2 = log
plotting_keys = [
    (keywords["id"], 0),
    (keywords["optimizer"], 1),
    (keywords["lr"], 2, [2e-4, 4e-1]),
    (keywords["initializer"], 1),
    (keywords["epochs"], 0),
    (keywords["ste"], 0),
    (keywords["stp"], 0),
    (keywords["p_mul"], 0),
    (keywords["nl"], 1),
    (keywords["activation"], 1),
    (keywords["dropout"], 0),
    (keywords["tl"], 0),
]


def parse_args():
    """ Wrapper around argumentparser """
    parser = ArgumentParser()
    parser.add_argument("folder", help="Fit folder")
    parser.add_argument(
        "-v",
        "--val_multiplier",
        help="Fraction to weight the validation loss with (test_multiplier = 1-val_multiplier)",
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "-if", "--include_failures", help="Flag to include failed runs in  the plots", action="store_true"
    )
    parser.add_argument(
        "-t",
        "--threshold",
        help="Value of the loss function from which to consider a run to have failed",
        type=float,
        default=1e3,
    )
    parser.add_argument("-f", "--filter", help="Add the filter key=value to the dataframe", nargs="+", default=())
    parser.add_argument(
        "-c", "--combine", help="If more than one replica folder is found, combine all trials", action="store_true"
    )
    # Autofiltering
    parser.add_argument(
        "--autofilter",
        help="Given a number of keys, perform an autofilter (removing combinations of elements with worse rewards",
        nargs="+",
    )
    args = parser.parse_args()
    return args


# Parse the different sections of the hyperoptimization
# this should talk with hyper_scan somehow?
def parse_optimizer(trial):
    """
    This function parses the parameters that affect the optimization

    optimizer
    learning_rate (if it exists)
    """
    dict_out = {}
    opt = trial["misc"]["space_vals"][keywords["optimizer"]]
    if isinstance(opt, dict):
        # If this is a dictionary then the optimizer contains extra
        # information (normaly the learning rate)
        name = opt[keywords["optimizer"]]
        lr = opt[keywords["lr"]]
    else:
        name = opt
        lr = None
    dict_out[keywords["optimizer"]] = name
    dict_out[keywords["lr"]] = lr
    return dict_out


def parse_stopping(trial):
    """
    This function parses the parameters that affect the stopping

    epochs
    stopping_patience
    pos_initial
    pos_multiplier
    """
    dict_out = {}
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
    return dict_out


def parse_architecture(trial):
    """
    This function parses the family of parameters which regards the architecture of the NN

    number_of_layers
    activation_per_layer
    nodes_per_layer
    l1, l2, l3, l4... max_layers
    layer_type
    dropout
    initializer
    """
    dict_out = {}
    nodes_per_layer = trial["misc"]["space_vals"][keywords["nodes"]]
    nl = len(nodes_per_layer) - 1
    activation_name = trial["misc"]["space_vals"][keywords["activation"]]
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

    ini = trial["misc"]["space_vals"][keywords["initializer"]]
    dropout_rate = trial["misc"]["space_vals"][keywords["dropout"]]
    dict_out[keywords["initializer"]] = ini
    dict_out[keywords["dropout"]] = dropout_rate
    return dict_out


def parse_statistics(trial):
    """
    Parse the statistical information of the trial

    validation loss
    testing loss
    status of the run
    """
    dict_out = {}
    validation_loss = trial["result"][keywords["vl"]]
    testing_loss = trial["result"][keywords["tl"]]
    # was this a ok run?
    ok = bool(trial["result"]["status"] == "ok")

    dict_out[keywords["good"]] = ok
    dict_out[keywords["vl"]] = validation_loss
    dict_out[keywords["tl"]] = testing_loss
    return dict_out


def parse_trial(trial):
    """
    Trials are very convoluted object, very branched inside
    The goal of this function is to separate said branching so we can create hierarchies
    """
    # Is this a true trial?
    if trial["state"] != 2:
        return None

    data_dict = {}
    # Parse all information into nicely looking dictionaries
    data_dict.update(parse_optimizer(trial))
    data_dict.update(parse_stopping(trial))
    data_dict.update(parse_architecture(trial))
    data_dict.update(parse_statistics(trial))
    return data_dict


def evaluate_trial(trial_dict, validation_multiplier, fail_threshold):
    """
    Read a trial dictionary and compute the true loss and decide whether the run passes or not
    """
    test_f = 1.0 - validation_multiplier
    val_loss = trial_dict[keywords["vl"]]
    test_loss = trial_dict[keywords["tl"]]
    loss = val_loss * validation_multiplier + test_loss * test_f

    if loss > fail_threshold or val_loss > fail_threshold or test_loss > fail_threshold:
        trial_dict["good"] = False
        # Set the loss an order of magnitude above the actual result so it shows obviously on the plots
        loss *= 10

    trial_dict["loss"] = loss


def generate_dictionary(
    replica_path, json_name="tries.json", starting_index=0, val_multiplier=0.5, fail_threshold=10.0
):
    """
        Reads a json file and returns a list of dictionaries

        # Arguments:
            - `replica_path`: folder in which the tries.json file can be found
            - `starting_index`: if the trials are to be added to an already existing
                                set, make sure the id has the correct index!
            - `val_multiplier`: validation multipler
            - `fail_threhsold`: threshold for the loss to consider a configuration as a failure
    """
    filename = "{0}/{1}".format(replica_path, json_name)

    # Open the file and reads it as a json
    with open(filename, "r") as jlist:
        input_trials = json.load(jlist)

    # Read all trials and create a list of dictionaries
    # which can be turn into a dataframe
    all_trials = []
    for tid, trial in enumerate(input_trials):
        index = starting_index + tid
        trial_dict = parse_trial(trial)
        if trial_dict is None:
            continue
        evaluate_trial(trial_dict, val_multiplier, fail_threshold)
        trial_dict[keywords["id"]] = index
        all_trials.append(trial_dict)
    return all_trials


def filter_by_string(filter_string):
    """
    Receives a data_dict (a parsed trial) and a filter string, returns True if the trial passes the filter

    filter string must have the format: key[]string
    where [] can be any of !=, =, >, <
    """

    def filter_function(data_dict):
        if filter_string is None:
            # No filter set, so everything passes
            return True

        match = regex_not_op.findall(filter_string)
        if len(match) < 2:
            raise ValueError("Filter str is not correct: {0}".format(filter_string))
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
            raise NotImplementedError("Filter string not valid, operator not recognized {0}".format(filter_string))

        # This I know it is not ok:
        if isinstance(val_check, str) and isinstance(filter_val, str):
            check_str = '"{0}"{1}"{2}"'
        else:
            check_str = "{0}{1}{2}"
        try:  # this whole thing is a bit naughty...
            return eval(check_str.format(val_check, operator, filter_val))
        except:
            return False

    return filter_function


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


def main():
    args = parse_args()

    # Prepare the filter(s)
    filter_functions = [filter_by_string(filter_me) for filter_me in args.filter]

    search_str = "{0}/nnfit/replica_*/tries.json".format(args.folder)
    all_json = glob.glob(search_str)
    starting_index = 0
    all_replicas = []
    for i, json_path in enumerate(all_json):
        # Look at the json and read all of them into a dictionary
        replica_path = os.path.dirname(json_path)
        dictionaries = generate_dictionary(
            replica_path,
            starting_index=starting_index,
            val_multiplier=args.val_multiplier,
            fail_threshold=args.threshold,
        )

        # Check if we are playing combinations, if we do keep reading json until we consume all of them
        if args.combine:
            starting_index += len(dictionaries)
            all_replicas += dictionaries
            # If this is not the last one, continue
            if (i + 1) == len(all_json):
                dictionaries = all_replicas
            else:
                continue

        # Now filter out the ones we don't want by passing every dictionary in the list
        # through all filters
        if filter_functions:
            valid_dictionaries = []
            for dictionary in dictionaries:
                if all(f(dictionary) for f in filter_functions):
                    valid_dictionaries.append(dictionary)
            dictionaries = valid_dictionaries

        # Now fill a pandas dataframe with the survivors of the filters
        dataframe_raw = pd.DataFrame(dictionaries)

        # If autofilter is active, apply it!
        if args.autofilter:
            name_keys = [keywords.get(i, i) for i in args.autofilter]
            dataframe_raw = autofilter_dataframe(dataframe_raw, name_keys)

        # By default we don't want failures
        if args.include_failures:
            dataframe = dataframe_raw
        else:
            dataframe = dataframe_raw[dataframe_raw["good"]]

        # Now select the best one
        best_idx = dataframe.loss.idxmin()
        best_trial_series = dataframe.ix[best_idx]
        # Make into a dataframe and transpose or the plotting code will complain
        best_trial = best_trial_series.to_frame().T

        # A script gotta plot what a script gotta plot
        fileout = "{0}/scan.pdf".format(replica_path)
        plot_scans(dataframe, best_trial, fileout)
        print("Plot saved at {0}".format(fileout))


if __name__ == "__main__":
    main()
