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
(and, of course the function in the fitting action that calls the minimization).
"""
import copy
import logging
from typing import Callable

import hyperopt
from hyperopt.pyll.base import scope
import numpy as np

from n3fit.backends import MetaLayer, MetaModel
from n3fit.hyper_optimization.filetrials import FileTrials
from n3fit.hyper_optimization.mongofiletrials import MongodRunner, MongoFileTrials

log = logging.getLogger(__name__)

# Hyperopt uses these strings for a passed and failed run
# it also has statuses "new", "running" and "suspended", but we don't use them
HYPEROPT_STATUSES = {True: "ok", False: "fail"}


HYPEROPT_SEED = 42


# These are just wrapper around some hyperopt's sampling expresions defined in here
# https://github.com/hyperopt/hyperopt/wiki/FMin#21-parameter-expressions
# with a bit of extra documentation for the ones that are not obvious
def hp_uniform(key, lower_end, higher_end):
    """Sample uniformly between lower_end and higher_end"""
    if lower_end is None or higher_end is None:
        return None
    return hyperopt.hp.uniform(key, lower_end, higher_end)


def hp_quniform(key, lower_end, higher_end, step_size=None, steps=None, make_int=False):
    """Like uniform but admits a step_size"""
    if lower_end is None or higher_end is None:
        return None
    if not step_size:
        step_size = lower_end
    if steps:
        step_size = (higher_end - lower_end) / steps

    ret = hyperopt.hp.quniform(key, lower_end, higher_end, step_size)
    if make_int:
        ret = scope.int(ret)
    return ret


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
    """Sample from the list or array ``choices``"""
    if len(choices) == 0:
        return None
    return hyperopt.hp.choice(key, choices)


# Wrapper for options that work the same
def optimizer_arg_wrapper(hp_key, option_dict):
    # Let's accept the learning rate as a fixed float
    if isinstance(option_dict, float):
        choice = option_dict
    else:
        # Get the min-max of the learing rate
        max_lr = option_dict["max"]
        min_lr = option_dict["min"]
        # Get the sampling type, if not given use uniform
        sampling = option_dict.get("sampling", "uniform")
        if sampling == "uniform":
            choice = hp_uniform(hp_key, min_lr, max_lr)
        elif sampling == "log":
            choice = hp_loguniform(hp_key, min_lr, max_lr)
    return choice


# Wrapper for the hyperscanning
def hyper_scan_wrapper(replica_path_set, model_trainer, hyperscanner, max_evals=1):
    """
    This function receives a ``ModelTrainer`` object as well as the definition of the
    hyperparameter scan (``hyperscanner``)
    and performs ``max_evals`` evaluations of the hyperparametrizable function of ``model_trainer``.

    A ``tries.json`` file will be saved in the ``replica_path_set`` folder with the information
    of all trials. An additional ``tries.pkl`` file will also be generated in the same folder
    that stores the previous states of `FileTrials`, this file can be used for restarting purposes.

    Parameters
    -----------
        replica_path_set: path
            folder where to create the ``tries.json`` and ``tries.pkl`` files
        model_trainer: :py:class:`n3fit.ModelTrainer.ModelTrainer`
            a ``ModelTrainer`` object with the ``hyperparametrizable`` method
        hyperscanner: :py:class:`n3fit.hyper_optimization.hyper_scan.HyperScanner`
            a ``HyperScanner`` object defining the scan
        max_evals: int
            Number of trials to run

    Returns
    -------
        dict
        parameters of the best trial as found by ``hyperopt``
    """
    # Tell the trainer we are doing hpyeropt
    model_trainer.set_hyperopt(True, keys=hyperscanner.hyper_keys)

    if hyperscanner.parallel_hyperopt:
        # start MongoDB database bu launching `mongod`
        hyperscanner.mongod_runner.ensure_database_dir_exists()
        mongod = hyperscanner.mongod_runner.start()

    if hyperscanner.restart_hyperopt:
        # For parallel hyperopt restarts, extract the database tar file
        if hyperscanner.parallel_hyperopt:
            tar_file_to_extract = f"{replica_path_set}/{hyperscanner.db_name}.tar.gz"
            log.info("Restarting hyperopt run using the MongoDB database %s", tar_file_to_extract)
            MongoFileTrials.extract_mongodb_database(tar_file_to_extract)
        else:
            # For sequential hyperopt restarts, reset the state of `FileTrials` saved in the pickle file
            pickle_file_to_load = f"{replica_path_set}/tries.pkl"
            log.info("Restarting hyperopt run using the pickle file %s", pickle_file_to_load)
            trials = FileTrials.from_pkl(pickle_file_to_load)

    # Generate the trials object
    if hyperscanner.parallel_hyperopt:
        # Instantiate `MongoFileTrials`
        # Mongo database should have already been initiated at this point
        trials = MongoFileTrials(
            replica_path_set,
            db_host=hyperscanner.db_host,
            db_port=hyperscanner.db_port,
            db_name=hyperscanner.db_name,
            num_workers=hyperscanner.num_mongo_workers,
            parameters=hyperscanner.as_dict(),
        )
    else:
        # Instantiate `FileTrials`
        trials = FileTrials(replica_path_set, parameters=hyperscanner.as_dict())

    # Initialize seed for hyperopt
    trials.rstate = np.random.default_rng(HYPEROPT_SEED)

    # Call to hyperopt.fmin
    fmin_args = dict(
        fn=model_trainer.hyperparametrizable,
        space=hyperscanner.as_dict(),
        algo=hyperopt.tpe.suggest,
        max_evals=max_evals,
        trials=trials,
        rstate=trials.rstate,
    )
    if hyperscanner.parallel_hyperopt:
        trials.start_mongo_workers()
        best = hyperopt.fmin(**fmin_args, show_progressbar=True, max_queue_len=trials.num_workers)
        trials.stop_mongo_workers()
        # stop mongod command and compress database
        hyperscanner.mongod_runner.stop(mongod)
        trials.compress_mongodb_database()
    else:
        best = hyperopt.fmin(**fmin_args, show_progressbar=False, trials_save_file=trials.pkl_file)
    return hyperscanner.space_eval(best)


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
        self._original_parameters = parameters
        self.parameter_keys = parameters.keys()
        self.parameters = copy.deepcopy(parameters)
        self.steps = steps

        # adding extra options for restarting
        restart_config = sampling_dict.get("restart")
        self.restart_hyperopt = True if restart_config else False

        # adding extra options for parallel execution
        parallel_config = sampling_dict.get("parallel")
        self.parallel_hyperopt = True if parallel_config else False

        # setting up MondoDB options
        if self.parallel_hyperopt:
            self.db_host = sampling_dict.get("db_host")
            self.db_port = sampling_dict.get("db_port")
            self.db_name = sampling_dict.get("db_name")
            self.num_mongo_workers = sampling_dict.get("num_mongo_workers")
            self.mongod_runner = MongodRunner(self.db_name, self.db_port)

        self.hyper_keys = set([])

        if "parameters" in sampling_dict:
            parameter_choices = hp_choice("parameters", sampling_dict["parameters"])
            # Drop the parameters dictionary
            self.parameters = {}
            self._update_param("parameters", parameter_choices)
            return

        # A hyperparameter scan will contain either a "parameters" or one of:
        stopping_dict = sampling_dict.get("stopping")
        optimizer_list = sampling_dict.get("optimizer")
        positivity_dict = sampling_dict.get("positivity")
        nn_dict = sampling_dict.get("architecture")

        if stopping_dict:
            self.stopping(
                min_epochs=stopping_dict.get("min_epochs"),
                max_epochs=stopping_dict.get("max_epochs"),
                min_patience=stopping_dict.get("min_patience"),
                max_patience=stopping_dict.get("max_patience"),
            )
        if optimizer_list:
            self.optimizer(optimizers=optimizer_list)
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
                output_size=parameters['nodes_per_layer'][-1],
            )

    def as_dict(self):
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

        if key not in self.parameter_keys and key != "parameters":
            raise ValueError(
                "Trying to update a parameter not declared in the `parameters` dictionary: {0} @ HyperScanner._update_param".format(
                    key
                )
            )

        self.hyper_keys.add(key)
        log.info("Adding key {0} with value {1}".format(key, sampler))

        self.parameters[key] = sampler

    def stopping(self, min_epochs=None, max_epochs=None, min_patience=None, max_patience=None):
        """
        Modifies the following entries of the `parameters` dictionary:
            - `epochs`
            - `stopping_patience`

        Takes `self.steps` between the min and maximum values given
        """
        epochs_key = "epochs"
        stopping_key = "stopping_patience"

        if min_epochs is not None and max_epochs is not None:
            epochs = hp_quniform(epochs_key, min_epochs, max_epochs, step_size=1, make_int=True)
            self._update_param(epochs_key, epochs)

        if min_patience is not None or max_patience is not None:
            if min_patience is None:
                min_patience = 0.0
            if max_patience is None:
                max_patience = 1.0

            stopping_patience = hp_quniform(
                stopping_key, min_patience, max_patience, steps=self.steps
            )

            self._update_param(stopping_key, stopping_patience)

    def optimizer(self, optimizers):
        """
        This function look at the optimizers implemented in MetaModel
        Since each optimizer can take different parameters, the input to this function, `optimizer`
        is a list of dictionaries, each defining the name of the optimizer (which needs to be
        implemented in `n3fit`) and the options to modify.

        The accepted options are:
            - learning_rate
            - clipnorm
        but for hyperopt it will look as a list of dictionaries
            [ { optimizer_name: optimizer_name, learning_rate: sampler },
              { optimizer_name: optimizer_name, learning_rate: sampler }, ...
            ]
        and will sample one from this list.

        Note that the keys within the dictionary (`optimizer_name` and `learning_rate`) should be named
        as the keys used by the compiler of the model as they are used as they come.
        """
        # Get all accepted optimizer to check against
        all_optimizers = MetaModel.accepted_optimizers
        # We will have a list of dictionaries to choose from
        choices = []

        opt_key = "optimizer"
        optname_key = "optimizer_name"
        lr_key = "learning_rate"
        clip_key = "clipnorm"

        for optimizer in optimizers:
            name = optimizer[optname_key]
            optimizer_dictionary = {optname_key: name}

            if name not in all_optimizers.keys():
                raise NotImplementedError(
                    f"HyperScanner: Optimizer {name} not implemented in MetaModel.py"
                )

            lr_dict = optimizer.get(lr_key)
            if lr_dict is not None:
                # Check whether this optimizer is implemented with a learning rate
                args = all_optimizers[name][1]
                if lr_key not in args.keys():
                    raise ValueError(f"Optimizer {name} does not accept {lr_key}")
                hp_key = f"{name}_{lr_key}"
                optimizer_dictionary[lr_key] = optimizer_arg_wrapper(hp_key, lr_dict)

            clip_dict = optimizer.get(clip_key)
            if clip_dict is not None:
                hp_key = f"{name}_{clip_key}"
                optimizer_dictionary[clip_key] = optimizer_arg_wrapper(hp_key, clip_dict)

            choices.append(optimizer_dictionary)

        # Make the list of options into a list sampler
        opt_val = hp_choice(opt_key, choices)

        # Tell the HyperScanner this key might contain a dictionary so we save the extra keys
        self._update_param(opt_key, opt_val)

    def positivity(
        self, min_multiplier=None, max_multiplier=None, min_initial=None, max_initial=None
    ):
        """
        Modifies the following entries of the `parameters` dictionary:
            - pos_multiplier
            - pos_initial
        Sampling between max and min is uniform for the multiplier and loguniform for the initial
        """
        mul_key = "multiplier"
        ini_key = "initial"
        params = {}

        if max_multiplier is not None:
            if min_multiplier is None:
                min_multiplier = 1.0  # I guess this is a sensible minimum

            mul_val = hp_uniform(mul_key, min_multiplier, max_multiplier)
            params[mul_key] = mul_val

        if min_initial is not None and max_initial is not None:
            ini_val = hp_loguniform(ini_key, min_initial, max_initial)
            params[ini_key] = ini_val

        self._update_param("positivity", params)

    def architecture(
        self,
        initializers=None,
        activations=None,
        max_drop=None,
        n_layers=None,
        min_units=15,
        max_units=25,
        layer_types=None,
        output_size=None,
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
                    units_label, min_units, max_units, step_size=1, make_int=True
                )
                units.append(units_sampler)
            # The number of nodes in the last layer are read from the runcard
            units.append(output_size)
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
                    "HyperScanner: Initializer {0} not implemented in MetaLayer.py".format(ini_name)
                )
            # For now we are going to use always all initializers and with default values
            ini_choices.append(ini_name)

        # Finally select the dropout rate, starting point always at 0

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

    def space_eval(self, trial):
        """Evaluate a trial using the original parameters dictionary"""
        return hyperopt.space_eval(self._original_parameters, trial)
