from argparse import ArgumentParser
import os
import re
import glob
import json
import logging
import numpy as np
import pandas as pd
from n3fit.hyper_optimization.hyper_algorithm import autofilter_dataframe
from types import SimpleNamespace

log = logging.getLogger(__name__)

regex_op = re.compile(r"[^\w^\.]+")
regex_not_op = re.compile(r"[\w\.]+")


# A mapping between the names the fields have in the json file
# and more user-friendly names
KEYWORDS = {
    "id": "iteration",
    "optimizer": "optimizer",
    "optimizer_name": "optimizer_name",
    "clipnorm": "clipnorm",
    "lr": "learning_rate",
    "initializer": "initializer",
    "dropout": "dropout",
    "nodes": "nodes_per_layer",
    "max_layers": 4,  # TODO: this should be dinamically choosen
    "nl": "number_of_layers",
    "activation": "activation_per_layer",
    "architecture": "layer_type",
    "epochs": "epochs",
    "stp": "stopping_patience",
    "ste": "stopping_epochs",
    "p_ini": "initial",
    "p_mul": "multiplier",
    "good": "good",
    "vl": "validation_loss",
    "tl": "loss",  # The testing loss has dissapeared, the loss corresponds to the k-folding loss
    "loss": "loss",
}


# Parse the different sections of the hyperoptimization
# this should talk with hyper_scan somehow?
def parse_optimizer(trial):
    """
    This function parses the parameters that affect the optimization

    optimizer
    learning_rate (if it exists)
    """
    dict_out = {}
    opt = trial["misc"]["space_vals"][KEYWORDS["optimizer"]]
    if isinstance(opt, dict):
        # If this is a dictionary then the optimizer contains extra
        # information (normaly the learning rate)
        name = opt[KEYWORDS["optimizer_name"]]
        lr = opt.get(KEYWORDS["lr"])
        clipnorm = opt.get(KEYWORDS["clipnorm"])
    else:
        name = opt
        lr = None
        clipnorm = None
    dict_out[KEYWORDS["optimizer"]] = name
    dict_out[KEYWORDS["lr"]] = lr
    dict_out[KEYWORDS["clipnorm"]] = clipnorm
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
    epochs = trial["misc"]["space_vals"][KEYWORDS["epochs"]]
    patience = trial["misc"]["space_vals"][KEYWORDS["stp"]]
    stop_ep = patience * epochs
    positivity_initial = trial["misc"]["space_vals"]["positivity"].get(KEYWORDS["p_ini"])
    positivity_multiplier = trial["misc"]["space_vals"]["positivity"].get(KEYWORDS["p_mul"])

    dict_out[KEYWORDS["epochs"]] = epochs
    dict_out[KEYWORDS["stp"]] = patience
    dict_out[KEYWORDS["ste"]] = stop_ep
    dict_out[KEYWORDS["p_ini"]] = positivity_initial
    dict_out[KEYWORDS["p_mul"]] = positivity_multiplier
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
    nodes_per_layer = trial["misc"]["space_vals"][KEYWORDS["nodes"]]
    nl = len(nodes_per_layer) - 1
    activation_name = trial["misc"]["space_vals"][KEYWORDS["activation"]]
    architecture = trial["misc"]["space_vals"][KEYWORDS["architecture"]]

    dict_out[KEYWORDS["nodes"]] = nodes_per_layer
    dict_out[KEYWORDS["nl"]] = nl
    dict_out[KEYWORDS["activation"]] = activation_name
    dict_out[KEYWORDS["architecture"]] = architecture

    # In principle we might be checking any number of layers,
    # so it is important to make sure that the maximum number of layers
    # (which will be used at plotting time) is correct
    if nl > KEYWORDS["max_layers"]:
        KEYWORDS["max_layers"] = nl

    for i in range(KEYWORDS["max_layers"]):
        dict_out["layer_{0}".format(i + 1)] = None

    for i, nodes in enumerate(nodes_per_layer[:-1]):
        dict_out["layer_{0}".format(i + 1)] = nodes

    ini = trial["misc"]["space_vals"][KEYWORDS["initializer"]]
    dropout_rate = trial["misc"]["space_vals"][KEYWORDS["dropout"]]
    dict_out[KEYWORDS["initializer"]] = ini
    dict_out[KEYWORDS["dropout"]] = dropout_rate
    return dict_out


def parse_statistics(trial):
    """
    Parse the statistical information of the trial

    validation loss
    testing loss
    status of the run
    """
    dict_out = {}
    results = trial["result"]
    validation_loss = results[KEYWORDS["vl"]]
    testing_loss = results[KEYWORDS["tl"]]
    # was this a ok run?
    ok = bool(results["status"] == "ok")

    dict_out[KEYWORDS["good"]] = ok
    dict_out[KEYWORDS["vl"]] = validation_loss
    dict_out[KEYWORDS["tl"]] = testing_loss

    # Kfolding information
    # average = results["kfold_meta"]["hyper_avg"]
    # std = results["kfold_meta"]["hyper_std"]
    # dict_out["avg"] = average
    # dict_out["std"] = std
    dict_out["hlosses"] = results["kfold_meta"]["hyper_losses"]
    dict_out["vlosses"] = results["kfold_meta"]["validation_losses"]
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
    val_loss = trial_dict[KEYWORDS["vl"]]
    test_loss = trial_dict[KEYWORDS["tl"]]
    loss = val_loss * validation_multiplier + test_loss * test_f

    if loss > fail_threshold or val_loss > fail_threshold or test_loss > fail_threshold:
        trial_dict["good"] = False
        # Set the loss an order of magnitude above the result so it shows obviously on the plots
        loss *= 10

    trial_dict["loss"] = loss


def generate_dictionary(
    replica_path,
    json_name="tries.json",
    starting_index=0,
    val_multiplier=0.5,
    fail_threshold=10.0,
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
        trial_dict[KEYWORDS["id"]] = index
        all_trials.append(trial_dict)
    return all_trials


def filter_by_string(filter_string):
    """
    Receives a data_dict (a parsed trial) and a filter string,
    returns True if the trial passes the filter

    filter string must have the format: key<operator>string
    where <operator> can be any of !=, =, >, <

    # Arguments:
        - `filter_string`: the expresion to evaluate

    # Returns:
        - `filter_function`: a function that takes a data_dict and
                             returns true if the condition in `filter_string` passes
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

        filter_key = KEYWORDS.get(filter_key, filter_key)

        val_check = data_dict[filter_key]
        if val_check is None:  # NaN means it does not apply
            return True

        operator = regex_op.findall(filter_string)[0]

        if operator == "=":
            operator = "=="
        operators = ["!=", "==", ">", "<"]
        if operator not in operators:
            raise NotImplementedError(
                "Filter string not valid, operator not recognized {0}".format(
                    filter_string
                )
            )

        # This I know it is not ok:
        if isinstance(val_check, str) and isinstance(filter_val, str):
            check_str = '"{0}"{1}"{2}"'
        else:
            check_str = "{0}{1}{2}"
        try:  # this whole thing is a bit naughty...
            return eval(check_str.format(val_check, operator, filter_val))
        except:  # if whatever happens within eval fails, just return False
            return False

    return filter_function


def hyperopt_dataframe(commandline_args):
    args = SimpleNamespace(**commandline_args)

    if args.debug:
        root_log = logging.getLogger()
        root_log.setLevel(logging.DEBUG)

    filter_functions = [filter_by_string(filter_me) for filter_me in args.filter]

    search_str = f"{args.hyperopt_folder}/nnfit/replica_*/tries.json"
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

        # Check if we are playing combinations,
        # if we do keep reading json until we consume all of them
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
            name_keys = [KEYWORDS.get(i, i) for i in args.autofilter]
            # First, for each key we are filtering in remove the worst
            # this will already remove a good chunk of trials
            for key in name_keys:
                dataframe_raw = autofilter_dataframe(
                    dataframe_raw, [key], n_to_combine=1, n_to_kill=1
                )
            how_many = len(name_keys)
            # Now remove the 2 worst for each combination from 2 to how_many
            for i in range(2, how_many + 1):
                dataframe_raw = autofilter_dataframe(
                    dataframe_raw, name_keys, n_to_combine=i, n_to_kill=2
                )

        # By default we don't want failures
        if args.include_failures:
            dataframe = dataframe_raw
        else:
            dataframe = dataframe_raw[dataframe_raw["good"]]

    # Now select the best one
    best_idx = dataframe.loss.idxmin()
    best_trial_series = dataframe.loc[best_idx]
    # Make into a dataframe and transpose or the plotting code will complain
    best_trial = best_trial_series.to_frame().T

    log.info("Best setup:")
    with pd.option_context("display.max_rows", None, "display.max_columns", None):
        log.info(best_trial)

    return dataframe, best_trial
