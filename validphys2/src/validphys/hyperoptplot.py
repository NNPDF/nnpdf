"""
    Module for the parsing and plotting of the results and output of
    previous hyperparameter scans
"""

# Within this file you can find the "more modern" vp-integrated hyperopt stuff
# and the older pre-vp hyperopt stuff, which can be considered deprecated but it is
# still used for the plotting script


import os
import re
import glob
import json
import logging
from types import SimpleNamespace
import numpy as np
import pandas as pd
from reportengine.figure import figure
from reportengine.table import table
import seaborn as sns
from matplotlib.figure import Figure
from validphys.hyper_algorithm import autofilter_dataframe

log = logging.getLogger(__name__)

regex_op = re.compile(r"[^\w^\.]+")
regex_not_op = re.compile(r"[\w\.]+")


class HyperoptTrial:
    """
    Hyperopt trial class.
    Makes the dictionary-like output of ``hyperopt`` into an object
    that can be easily managed

    Parameters
    ----------
        trial_dict: dict
            one single result (a dictionary) from a ``tries.json`` file
        base_params: dict
            Base parameters of the runcard which can be used to complete the hyperparameter
            dictionary when not all parameters were scanned
        minimum_losses: int
            Minimum number of losses to be found in the trial for it to be considered succesful
        linked_trials: list
            List of trials coming from the same file as this trial
    """

    def __init__(self, trial_dict, base_params=None, minimum_losses=1, linked_trials=None):
        self._trial_dict = trial_dict
        self._minimum_losses = minimum_losses
        self._original_params = base_params
        self._reward = None
        self._weighted_reward = None
        self._linked_trials = linked_trials if linked_trials is not None else []

    @property
    def weighted_reward(self):
        """Return the reward weighted to the mean value of the linked trials"""
        if self._weighted_reward is None:
            mean_reward = np.mean([i.reward for i in self._linked_trials if i.reward])
            if self.reward:
                self._weighted_reward = self.reward / mean_reward
            else:
                self._weighted_reward = False
        return self._weighted_reward

    @property
    def reward(self):
        """Return and cache the reward value"""
        if self._reward is None:
            self._reward = self._evaluate()
        return self._reward

    @property
    def loss(self):
        """Return the loss of the hyperopt dict"""
        return self._trial_dict["result"]["loss"]

    @property
    def params(self):
        """Parameters for the fit"""
        trial_results = self._trial_dict["misc"]["space_vals"]
        if self._original_params:
            hyperparameters = dict(self._original_params, **trial_results)
        else:
            hyperparameters = trial_results
        # Ensure that no hyperparameters has the wrong shape/form
        # 1. For n3fit, activation_per_layer must be a list
        apl = hyperparameters.get("activation_per_layer")
        if apl is not None and isinstance(apl, str) and "nodes_per_layer" in hyperparameters:
            # If apl is a string, it can only bring information if there is `nodes_per_layer`
            apl = [apl] * (len(hyperparameters["nodes_per_layer"]) - 1) + ["linear"]
            hyperparameters["activation_per_layer"] = apl
        # 2. If there _was_ originally a reward included in this parameter (because it comes from a
        #   previous hyperoptimization) mark it as such
        reward = hyperparameters.pop("reward", None)
        if reward is not None:
            hyperparameters["old_reward"] = reward
        return hyperparameters

    # Slightly fake a dictionary behaviour with the params attribute
    def __getitem__(self, item):
        return self.params[item]

    def get(self, item, default=None):
        return self.params.get(item, default)

    #####

    def __str__(self):
        strs = ["Parameters:"] + [f"  {i}: {k}" for i, k in self.params.items()]
        strs.append(f"Reward: {self.reward}")
        str_out = "\n".join(strs)
        return str_out

    def _evaluate(self):
        """Evaluate the reward function for a given trial
        Reward = 1.0/loss
        """
        result = self._trial_dict["result"]
        if result.get("status") != "ok":
            return False
        # If we don't have enough validation losses, fail
        val_loss = result["kfold_meta"].get("validation_losses", [])
        if self._minimum_losses and len(val_loss) < self._minimum_losses:
            return False
        # This is equivalent to the case above
        if result["loss"] == 0.0:
            return False
        return 1.0 / result["loss"]

    def __gt__(self, another_trial):
        """Return true if the current trial is preferred
        when compared to the target"""
        if not another_trial.reward:
            return True
        if not self.reward:
            return False
        return self.reward > another_trial.reward

    def __lt__(self, another_trial):
        """Return true if the target trial is preferred
        when compared to the current (self) one"""
        return not self > another_trial

    def link_trials(self, list_of_trials):
        """Link a list of trials to this trial"""
        self._linked_trials = list_of_trials
        # Reset the weighted reward
        self._weighted_reward = None


#########

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


# 0 = normal scatter plot, 1 = violin, 2 = log
plotting_styles = {
    "iteration": 0,
    "optimizer": 1,
    "learning_rate": 2,
    "initializer": 1,
    "epochs": 0,
    "stopping_epochs": 0,
    "stopping_patience": 0,
    "multiplier": 0,
    "number_of_layers": 1,
    "activation_per_layer": 1,
    "dropout": 0,
    "clipnorm": 0,
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


def evaluate_trial(trial_dict, validation_multiplier, fail_threshold, loss_target):
    """
    Read a trial dictionary and compute the true loss and decide whether the run passes or not
    """
    test_f = 1.0 - validation_multiplier
    val_loss = float(trial_dict[KEYWORDS["vl"]])
    if loss_target == "average":
        test_loss = np.array(trial_dict["hlosses"]).mean()
    elif loss_target == "best_worst":
        test_loss = np.array(trial_dict["hlosses"]).max()
    elif loss_target == "std":
        test_loss = np.array(trial_dict["hlosses"]).std()
    loss = val_loss * validation_multiplier + test_loss * test_f

    if (
        loss > fail_threshold
        or val_loss > fail_threshold
        or test_loss > fail_threshold
        or np.isnan(loss)
    ):
        trial_dict["good"] = False
        # Set the loss an order of magnitude above the result so it shows obviously on the plots
        loss *= 10

    trial_dict["loss"] = loss


def generate_dictionary(
    replica_path,
    loss_target,
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
        evaluate_trial(trial_dict, val_multiplier, fail_threshold, loss_target)
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
                "Filter string not valid, operator not recognized {0}".format(filter_string)
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
    """
    Loads the data generated by running hyperopt and stored in json files into a dataframe, and
    then filters the data according to the selection criteria provided by the command line
    arguments. It then returns both the entire dataframe as well as a dataframe object with the
    hyperopt parametesr of the best setup.
    """
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
            args.loss_target,
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


@table
def best_setup(hyperopt_dataframe, hyperscan_config, commandline_args):
    """
    Generates a clean table with information on the hyperparameter settings of the best setup.
    """
    _, best_trial = hyperopt_dataframe
    best_idx = best_trial.index[0]
    best_trial = best_trial.rename(index={best_idx: "parameter settings"})
    best_trial = best_trial[
        [
            "optimizer",
            "learning_rate",
            "clipnorm",
            "epochs",
            "stopping_patience",
            "initial",
            "multiplier",
            "nodes_per_layer",
            "activation_per_layer",
            "initializer",
            "dropout",
            "loss",
        ]
    ]
    best_trial.insert(11, "loss type", commandline_args["loss_target"])
    best_trial = best_trial.T
    return best_trial


@table
def hyperopt_table(hyperopt_dataframe):
    """
    Generates a table containing complete information on all the tested setups that passed the
    filters set in the commandline arguments.
    """
    dataframe, _ = hyperopt_dataframe
    dataframe.sort_values(by=["loss"], inplace=True)
    return dataframe


@figure
def plot_iterations(hyperopt_dataframe):
    """
    Generates a scatter plot of the loss as a function of the iteration index.
    """
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "iteration")
    return fig


@figure
def plot_optimizers(hyperopt_dataframe):
    """
    Generates a violin plot of the loss per optimizer.
    """
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "optimizer")
    return fig


@figure
def plot_clipnorm(hyperopt_dataframe, optimizer_name):
    """
    Generates a scatter plot of the loss as a function of the clipnorm for a given optimizer.
    """
    dataframe, best_trial = hyperopt_dataframe
    filtered_dataframe = dataframe[dataframe.optimizer == optimizer_name]
    best_filtered_idx = filtered_dataframe.loss.idxmin()
    best_idx = best_trial.iteration.iloc[0]
    if best_filtered_idx == best_idx:
        include_best = True
    else:
        include_best = False
    fig = plot_scans(filtered_dataframe, best_trial, "clipnorm", include_best=include_best)
    return fig


@figure
def plot_learning_rate(hyperopt_dataframe, optimizer_name):
    """
    Generates a scatter plot of the loss as a function of the learning rate for a given optimizer.
    """
    dataframe, best_trial = hyperopt_dataframe
    filtered_dataframe = dataframe[dataframe.optimizer == optimizer_name]
    best_filtered_idx = filtered_dataframe.loss.idxmin()
    best_idx = best_trial.iteration.iloc[0]
    if best_filtered_idx == best_idx:
        include_best = True
    else:
        include_best = False
    fig = plot_scans(filtered_dataframe, best_trial, "learning_rate", include_best=include_best)
    return fig


@figure
def plot_initializer(hyperopt_dataframe):
    """
    Generates a violin plot of the loss per initializer.
    """
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "initializer")
    return fig


@figure
def plot_epochs(hyperopt_dataframe):
    """
    Generates a scatter plot of the loss as a function the number of epochs.
    """
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "epochs")
    return fig


@figure
def plot_number_of_layers(hyperopt_dataframe):
    """
    Generates a violin plot of the loss as a function of the number of layers of the model.
    """
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "number_of_layers")
    return fig


@figure
def plot_activation_per_layer(hyperopt_dataframe):
    """
    Generates a violin plot of the loss per activation function.
    """
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "activation_per_layer")
    return fig


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


def plot_scans(df, best_df, plotting_parameter, include_best=True):
    """
    This function performs the plotting and is called by the `plot_` functions in this file.
    """
    figs = Figure()
    ax = figs.subplots()

    # Set the quantity we will be plotting in the y axis
    loss_k = "loss"

    key = plotting_parameter
    mode = plotting_styles[plotting_parameter]

    if mode in (0,2):  # normal scatter plot
        ax = sns.scatterplot(x=key, y=loss_k, data=df, ax=ax)
        best_x = best_df.get(key)
        if mode == 2:
            ax.set_xscale("log")
    elif mode == 1:
        sample = best_df.get(key).tolist()[0]
        if isinstance(sample, list):
            # activation_per_layer is tricky as it can be a list (with the last layer linear)
            # and can change size, the legacy way of plotting it was to take just the first function
            # For that we'll modify the dataframe that we pass down
            original_column = df[key]
            original_best = best_df[key]
            key += "_0"
            if key not in df:
                df[key] = original_column.apply(lambda x: x[0])
                best_df[key] = original_best.apply(lambda x: x[0])
        ordering_true, best_x = order_axis(df, best_df, key=key)
        ax = sns.violinplot(
            x=key,
            y=loss_k,
            data=df,
            ax=ax,
            palette="Set2",
            cut=0.0,
            order=ordering_true,
        )
        ax = sns.stripplot(
            x=key,
            y=loss_k,
            data=df,
            ax=ax,
            color="gray",
            order=ordering_true,
            alpha=0.4,
        )


    # Finally plot the "best" one, which will be first
    if include_best:
        ax = sns.scatterplot(x=best_x, y=best_df.get(loss_k), ax=ax, color="orange", marker="s")
    ax.set_ylabel("Loss")
    ax.set_xlabel(key)

    return figs
