"""
The HyperScanner class is basically a dictionary containing all parameters,
the functions defined as hp_ (from hyperspace)

The goal of this module is to read all parameters in the `hyperopt` section of the runcard
and modify the parameter dictionary so that now it is filled with the `hyperop` sampler objects

The idea behind the wrappers if that if you ever want to use
another hyperoptimization library, assuming that it also takes just
    - a function
    - a dictionary of spaces of parameters
you can do so by simply modifying the wrappers to point somewhere else
(and, of course the function in the fitting action that calls the miniimization).
"""
# TODO: make this part of Report Engine
# in principle Report Engine would be reading everything from the runcard and then providing
# basically a list of strings or similar to ModelTrainer, which would in turn make them into
# the keras/tensorflow/whatever objects it needs.
# If the --hpyeropt flag is active, then report engine should return a sampler from the dict
# basically returning the wrapper functions defined below
# (but instead of simply calling the hyperopt ones having done some checks on them)
# even more, depending on the syntax of the fitting part of the runcard,
# it should do one thing or another but the design of this needs someone who is actually good at UX

import copy
import hyperopt
import numpy as np
from n3fit.backends import MetaModel, MetaLayer
import n3fit.hyper_optimization.filetrials as filetrials
import logging

log = logging.getLogger(__name__)

# These are just wrapper around some hyperopt's sampling expresions defined in here
# https://github.com/hyperopt/hyperopt/wiki/FMin#21-parameter-expressions
# with a bit of extra documentation for the ones that are not obvious
def hp_uniform(key, lower_end, higher_end):
    """ Sample uniformly between lower_end and higher_end """
    if lower_end is None or higher_end is None:
        return None
    return hyperopt.hp.uniform(key, lower_end, higher_end)


def hp_quniform(key, lower_end, higher_end, step_size=None, steps=None):
    """ Like uniform but admits a step_size """
    if lower_end is None or higher_end is None:
        return None
    if not step_size:
        step_size = lower_end
    if steps:
        step_size = (higher_end - lower_end) / steps
    return hyperopt.hp.quniform(key, lower_end, higher_end, step_size)


def hp_loguniform(key, lower_end, higher_end):
    """
    Sample from lower_end to higher_end lograithmically.
    Note that it is different from numpy's logspace in that it takes the
    lower and higher boundaries, not the value of the exponent
    """
    if lower_end is None or higher_end is None:
        return None
    log_lower_end = np.log(lower_end)
    log_higher_end = np.log(higher_end)
    return hyperopt.hp.loguniform(key, log_lower_end, log_higher_end)


def hp_choice(key, choices):
    """ Sample from the list `choices` """
    if not choices:
        return None
    return hyperopt.hp.choice(key, choices)


# Wrapper for the hyperscanning
def hyper_scan_wrapper(
    replica_path_set, model_trainer, parameters, hyperscan_dict, max_evals=1
):
    """
    This function receives a `ModelTrainer` object as well as a dictionary of parameters (`parameters`) and a dictionary defining the
    space to search (`hyperscan_dict`) and performs `max_evals` evaluations of the function trying to find the best model.

    The return value of this function is a dictionary containing a dictionary similar to `parameters` with the best model found.
    It will also write a json file with the information of all separate trials inside the replica_path_set directory.

    # Arguments:
        - `replica_path_set`: folder where to create json file with info about all trials
        - `model_trainer`: a ModelTrainer object
        - `parameters`: a dictionary of parameters to pass to hyperparametizable
        - `hyperscan_dict`: dictionary definining the sampling
        - `max_evals`: how many runs of hyperopt to run
    """
    # Create a HyperScanner object
    the_scanner = HyperScanner(parameters, hyperscan_dict)
    # Tell the trainer we are doing hpyeropt
    model_trainer.set_hyperopt(
        True, keys=the_scanner.hyper_keys, status_ok=hyperopt.STATUS_OK
    )
    # Generate the trials object
    trials = filetrials.FileTrials(replica_path_set, parameters=the_scanner.dict())

    # Perform the scan
    best = hyperopt.fmin(
        fn=model_trainer.hyperparametizable,
        space=the_scanner.dict(),
        algo=hyperopt.tpe.suggest,
        max_evals=max_evals,
        show_progressbar=False,
        trials=trials,
    )
    true_best = hyperopt.space_eval(parameters, best)
    return true_best


###


class ActivationStr:
    """
    Upon call this class returns an array where the activation function
    `fun_name` repeated as many times as hidden layers we have

    # Arguments:
        - `fun_name`: name of the activation function
    """

    def __init__(self, fun_name):
        self.function_name = fun_name

    def __str__(self):
        return self.function_name

    def __call__(self, n_of_layers):
        acts = [self.function_name] * (n_of_layers - 1)
        acts.append("linear")
        return acts


class HyperScanner:
    """
    The HyperScanner generates a dictionary of parameters for scanning
    It takes cares of known correlation between parameters by tying them together
    It also provides methods for updating the parameter dictionaries after using hyperopt

    It takes as inpujt the dictionaries defining the NN/fit and the hyperparameter scan
    from the NNPDF runcard and substitutes in `parameters` samplers according to the
    `hyper_scan` dictionary.


    # Arguments:
        - `parameters`: the `fitting[parameters]` dictionary of the NNPDF runcard
        - `sampling_dict`: the `hyperscan` dictionary of the NNPDF runcard defining
                           the search space of the scan
        - `steps`: when taking discrete steps between two parameters, number of steps
                   to take

    # Parameters accepted by `sampling_dict`:
        - `stopping`:
                - min_epochs, max_epochs
                - min_patience, max_patience
    """

    def __init__(self, parameters, sampling_dict, steps=5):
        self.parameter_keys = parameters.keys()
        self.parameters = copy.deepcopy(parameters)
        self.steps = steps

        self.hyper_keys = set([])

        stopping_dict = sampling_dict.get("stopping")
        optimizer_dict = sampling_dict.get("optimizer")
        positivity_dict = sampling_dict.get("positivity")
        nn_dict = sampling_dict.get("architecture")

        if stopping_dict:
            self.stopping(
                min_epochs=stopping_dict.get("min_epochs"),
                max_epochs=stopping_dict.get("max_epochs"),
                min_patience=stopping_dict.get("min_patience"),
                max_patience=stopping_dict.get("max_patience"),
            )
        if optimizer_dict:
            self.optimizer(
                names=optimizer_dict.get("names"),
                min_lr=optimizer_dict.get("min_lr"),
                max_lr=optimizer_dict.get("max_lr"),
            )
        if positivity_dict:
            self.positivity(
                min_multiplier=positivity_dict.get("min_multiplier"),
                max_multiplier=positivity_dict.get("max_multiplier"),
                min_initial=positivity_dict.get("min_initial"),
                max_initial=positivity_dict.get("max_initial"),
            )
        if nn_dict:
            self.architecture(
                initializers=nn_dict.get("initializers"),
                activations=nn_dict.get("activations"),
                max_drop=nn_dict.get("max_drop"),
                n_layers=nn_dict.get("n_layers"),
                min_units=nn_dict.get("min_units"),
                max_units=nn_dict.get("max_units"),
                layer_types=nn_dict.get("layer_types"),
            )

    def dict(self):
        return self.parameters

    def _update_param(self, key, sampler):
        """
        Checks whether the key exists in the parameter dictionary and
        updates the dictionary with the given sampler

        # Arguments:
            - `key`: key to update
            - `sampler`: sampler which will be used instead of the original value
        """
        if key is None or sampler is None:
            return

        if key not in self.parameter_keys:
            raise ValueError(
                "Trying to update a parameter not declared in the `parameters` dictionary: {0} @ HyperScanner._update_param".format(
                    key
                )
            )

        self.hyper_keys.add(key)
        log.info("Adding key {0} with value {1}".format(key, sampler))

        self.parameters[key] = sampler

    def stopping(
        self, min_epochs=5e3, max_epochs=30e3, min_patience=0.10, max_patience=0.3
    ):
        """
        Modifies the following entries of the `parameters` dictionary:
            - `epochs`
            - `stopping_patience`

        Takes `self.steps` between the min and maximum values given
        """
        epochs_key = "epochs"
        stopping_key = "stopping_patience"

        # Generate the samplers
        epochs = hp_quniform(epochs_key, min_epochs, max_epochs, steps=self.steps)
        stopping_patience = hp_quniform(
            stopping_key, min_patience, max_patience, steps=self.steps
        )

        # Update the parameters ditionary
        self._update_param(epochs_key, epochs)
        self._update_param(stopping_key, stopping_patience)

    def optimizer(self, names=None, min_lr=0.0005, max_lr=0.5):
        """
        This function look at the optimizers implemented in MetaModel and adds the learning rate to (only)
        those who use it. The special keyword "ALL" will make it use all optimizers implemented
            - `optimizer`
            - `learning rate`

        Since the learning rate is a parameter that depends on the optimizer, the final result is a
        recursive dictionary such that even though the ModelTrainer will receive
            {optimizer: sampler(optimizer), learning rate: sampler(lr)}
        for hyperopt it will look as
            {optimizer: [ (optimizer1, sampler(lr)), (optimzier2, sampler(lr)), (optimizer3, )]
        """
        if names is None:  # Nothing to do here
            return

        opt_key = "optimizer"
        lr_key = "learning_rate"

        # Check which optimizers are we using
        optimizer_dict = MetaModel.accepted_optimizers
        if names == "ALL":
            names = optimizer_dict.keys()

        # Set a logarithmic sampling for the learning rate
        lr_choice = hp_uniform(lr_key, min_lr, max_lr)

        choices = []
        for opt_name in names:
            if opt_name not in optimizer_dict.keys():
                raise NotImplementedError(
                    "HyperScanner: Optimizer {0} not implemented in MetaModel.py".format(
                        opt_name
                    )
                )

            # Check whether this optimizer takes the learning rate
            # if it does the choice should be a dictionary with both the optimizer and the learning rate
            # otherwise just the optimizer name (dictionary not needed)
            args = optimizer_dict[opt_name][1]
            if "lr" in args.keys():
                choices.append({opt_key: opt_name, lr_key: lr_choice})
            else:
                choices.append(opt_name)

        # Make the list of options into a list sampler
        opt_val = hp_choice(opt_key, choices)

        # Tell the HyperScanner this key might contain a dictionary so we save the extra keys
        self._update_param(opt_key, opt_val)

    def positivity(
        self, min_multiplier=1.01, max_multiplier=1.3, min_initial=0.5, max_initial=100
    ):
        """
        Modifies the following entries of the `parameters` dictionary:
            - pos_multiplier
            - pos_initial
        Sampling between maximum and minimum is uniform for the multiplier and loguniform for the initial
        """
        mul_key = "pos_multiplier"
        ini_key = "pos_initial"

        # Create the samplers
        mul_val = hp_uniform(mul_key, min_multiplier, max_multiplier)
        ini_val = hp_loguniform(ini_key, min_initial, max_initial)

        # Update the dictionaries
        self._update_param(mul_key, mul_val)
        self._update_param(ini_key, ini_val)

    def architecture(
        self,
        initializers=None,
        activations=None,
        max_drop=None,
        n_layers=None,
        min_units=15,
        max_units=25,
        layer_types=None,
    ):
        """
        Modifies the following entries of the `parameters` dictionary:
            - `initializer`
            - `dropout`
            - `nodes_per_layer`
            - `activation_per_layer`
            - `layer_type`
        """
        if activations is None:
            activations = []
        if initializers is None:
            initializers = []
        if layer_types is None:
            layer_types = []
        if n_layers is None:
            n_layers = []
        else:
            if min_units is None or max_units is None:
                raise ValueError(
                    "A max/min number of units must always be defined if the number of layers is to be sampled"
                    "i.e., make sure you add the keywords 'min_units' and 'max_units' to the 'architecutre' dict"
                )

        activation_key = "activation_per_layer"
        nodes_key = "nodes_per_layer"
        ini_key = "initializer"

        # Generate all possible activation choices
        activation_choices = []
        for afun in activations:
            activation_str = ActivationStr(afun)
            activation_choices.append(activation_str)

        # this is strongly coupled with the total number of layers
        # so we will generate a list of layers to choose from
        # where each layer will be defined by an uniform sampler (the number of nodes)
        nodes_choices = []
        for n in n_layers:
            units = []
            for i in range(n):
                units_label = "nl{0}:-{1}/{0}".format(n, i)
                units_sampler = hp_quniform(
                    units_label, min_units, max_units, steps=self.steps
                )
                units.append(units_sampler)
            # The last layer will always have 8 nodes
            units.append(8)
            nodes_choices.append(units)

        # For the initializer we need to check for the ones implemented in MetaLayer
        imp_inits = MetaLayer.initializers
        imp_init_names = imp_inits.keys()
        if initializers == "ALL":
            initializers = imp_init_names

        ini_choices = []
        for ini_name in initializers:
            if ini_name not in imp_init_names:
                raise NotImplementedError(
                    "HyperScanner: Initializer {0} not implemented in MetaLayer.py".format(
                        ini_name
                    )
                )
            # For now we are going to use always all initializers and with default values
            ini_choices.append(ini_name)


        # Create the samplers
        act_functions = hp_choice(activation_key, activation_choices)
        nodes = hp_choice(nodes_key, nodes_choices)
        ini_choice = hp_choice(ini_key, ini_choices)


        # Finally select the layer types (not very well tested for now)
        layer_key = "layer_type"
        layer_choices = hp_choice(layer_key, layer_types)
        self._update_param(layer_key, layer_choices)

        # And update the dictionary
        self._update_param(activation_key, act_functions)
        self._update_param(nodes_key, nodes)
        self._update_param(ini_key, ini_choice)

        if max_drop is not None:
            # Finally select the dropout rate, starting point always at 0
            drop_key = "dropout"
            drop_val = hp_quniform(drop_key, 0.0, max_drop, steps=self.steps)
            self._update_param(drop_key, drop_val)
