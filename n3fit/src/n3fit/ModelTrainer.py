"""
    The ModelTrainer class is the true driver around the n3fit code

    This class is initialized with all information about the NN, inputs and outputs.
    The construction of the NN and the fitting is performed at the same time when the
    hyperparametrizable method of the function is called.

    This allows to use hyperscanning libraries, that need to change the parameters of the network
    between iterations while at the same time keeping the amount of redundant calls to a minimum
"""
import logging
import numpy as np
import n3fit.model_gen as model_gen
from n3fit.backends import MetaModel, clear_backend_state, operations, callbacks
from n3fit.stopping import Stopping
import n3fit.hyper_optimization.penalties
import n3fit.hyper_optimization.rewards

log = logging.getLogger(__name__)

# Threshold defaults
# Any partition with a chi2 over the threshold will discard its hyperparameters
HYPER_THRESHOLD = 50.0
CHI2_THRESHOLD = 10.0
# Each how many epochs do we increase the positivitiy Lagrange Multiplier
PUSH_POSITIVITY_EACH = 100

# Each how many epochs do we increase the integrability Lagrange Multiplier
PUSH_INTEGRABILITY_EACH = 100


def _assign_data_to_model(model, data_dict, fold_k=0):
    """
    Reads the data dictionary (``data_dict``) and assings the target data to the model
    It returns a dictionary containing:
    {
        'model': the backend.MetaModel to be trained
        'target_ndata': an array of target output
        'ndata': the number of datapoints
        'losses': the list of loss functions of the model
    }

    If kfolding is active applies the (``fold_k``) fold to the target data.
    in this case ndata is a count of the non_zero entries of the fold

    Note: this function is transitional. Eventually a validphys action should
    provide an experiment object with a .target_data(fold_indx) method which
    should return the necessary information:
        - number of datapoints
        - target data with the right entries set to 0*
    *or masked away if it is able to also return a list of loss functions that will mask away
    the corresponding entries of the prediction


    Parameters
    ----------
        model: backend.MetaModel
            model to be added to the dictionary
        data_dict: dict
            dictionary containing: {
                'expdata' : list of experimental data which the model will target,
                'folds' : a list (size=expdata) of lists (size=kfolds) with the folding masks)
                'losses': a list of loss functions for the model
                }
        fold_k: int
            when kfolding, index of the fold, so that for every experiment we apply the
            folds[index_experiment][mask]

    Returns
    -------
        ret: dict
            dictionary containing the model and its associated ndata and loss
    """
    # Array with all data
    all_data = data_dict["expdata"]
    # Each element of this list correspond to the set of folds for one experiment
    all_folds = data_dict["folds"]
    n_exps = len(all_folds)
    # Now set to 0 the data folded away
    active_data = []
    ndata = 0
    for exp_data, exp_fold in zip(all_data, all_folds):
        if exp_fold:
            mask = exp_fold[fold_k]
            active_data.append(exp_data * mask)
            ndata += np.count_nonzero(mask)
        else:
            active_data.append(exp_data)
            ndata += exp_data.size
    # There might be special outputs (like positivitiy) that is not
    # affected by the folding.
    # They don't count for the chi2 (which is only for reporting)
    active_data += all_data[n_exps:]
    ret = {
        "model": model,
        "target_data": active_data,
        "ndata": ndata,
        "losses": data_dict["losses"],
    }
    return ret


def _model_compilation(models, optimizer_params):
    """
        Compiles all models

    Parameters
    ----------
        models: list(dict)
            A ditionary defining the model
        optimizer_params: dict
            Optimizer parameters to be passes to the compile method of the model
    """
    for _, model_dict in models.items():
        model = model_dict["model"]
        target = model_dict["target_data"]
        losses = model_dict["losses"]
        model.compile(loss=losses, target_output=target, **optimizer_params)


def _pdf_injection(pdf_layers, observables, datasets_out=None):
    """
    Takes as input a list of output layers and returns a corresponding list
    where all output layers call the pdf layer at self.pdf_layer
    """
    return [f(x, datasets_out=datasets_out) for f, x in zip(observables, pdf_layers)]


def _LM_initial_and_multiplier(input_initial, input_multiplier, max_lambda, steps):
    """
    If any of input_initial or input_multiplier is None this function computes
    the missing values taking as input the maximum lambda multiplier and the number of steps needed
    to reach the maximum number of epochs
    """
    initial = input_initial
    multiplier = input_multiplier
    # If the multiplier is None, compute it from known values
    if multiplier is None:
        # If the initial value is also None, set it to one
        if initial is None:
            initial = 1.0
        multiplier = pow(max_lambda / initial, 1 / max(steps, 1))
    elif initial is None:
        # Select the necessary initial value to get to max_lambda after all steps
        initial = max_lambda / pow(multiplier, steps)
    return initial, multiplier


class ModelTrainer:
    """
    ModelTrainer Class:

    Wrapper around the fitting code and the generation of the Neural Network

    When the "hyperparametrizable"* function is called with a dictionary of parameters,
    it generates a NN and subsequentially performs a fit.

    The motivation behind this class is minimising the amount
    of redundant calls of each hyperopt run, in particular this allows to completely reset
    the NN at the beginning of each iteration reusing some of the previous work.

    *called in this way because it accept a dictionary of hyper-parameters
    which defines the Neural Network
    """

    def __init__(
        self,
        exp_info,
        pos_info,
        integ_info,
        flavinfo,
        fitbasis,
        nnseed,
        pass_status="ok",
        failed_status="fail",
        debug=False,
        save_weights_each=False,
        kfold_parameters=None,
        max_cores=None,
        model_file=None, 
    ):
        """
        Parameters
        ----------
            exp_info: list of dictionaries containing experiments
            pos_info: list of dictionaries containing positivity sets
            integ_info: list of dictionaries containing integrability sets
            flavinfo: the object returned by fitting['basis']
            nnseed: the seed used to initialise the Neural Network, will be passed to model_gen
            pass_status: flag to signal a good run
            failed_status: flag to signal a bad run
            debug: flag to activate some debug options
            save_weights_each: if set, save the state of the fit
                                    every ``save_weights_each`` epochs
            model_file: str whether to save the models
        """

        # Save all input information
        self.exp_info = exp_info
        self.pos_info = pos_info
        self.integ_info = integ_info
        if self.integ_info is not None:
            self.all_info = exp_info + pos_info + integ_info
        else:
            self.all_info = exp_info + pos_info
        self.flavinfo = flavinfo
        self.fitbasis = fitbasis
        self.NNseed = nnseed
        self.pass_status = pass_status
        self.failed_status = failed_status
        self.debug = debug
        self.save_weights_each = save_weights_each
        self.all_datasets = []

        # Initialise internal variables which define behaviour
        if debug:
            self.max_cores = 1
        else:
            self.max_cores = max_cores
        self.model_file = model_file
        self.print_summary = True
        self.mode_hyperopt = False
        self.impose_sumrule = True
        self.hyperkeys = None
        if kfold_parameters is None:
            self.kpartitions = [None]
            self.hyper_threshold = None
        else:
            self.kpartitions = kfold_parameters["partitions"]
            self.hyper_threshold = kfold_parameters.get("threshold", HYPER_THRESHOLD)
            # if there are penalties enabled, set them up
            penalties = kfold_parameters.get("penalties", [])
            self.hyper_penalties = []
            for penalty in penalties:
                pen_fun = getattr(n3fit.hyper_optimization.penalties, penalty)
                self.hyper_penalties.append(pen_fun)
                log.info("Adding penalty: %s", penalty)
            # Check what is the hyperoptimization target function
            hyper_loss = kfold_parameters.get("target", None)
            if hyper_loss is None:
                hyper_loss = "average"
                log.warning("No minimization target selected, defaulting to '%s'", hyper_loss)
            log.info("Using '%s' as the target for hyperoptimization", hyper_loss)
            self.hyper_loss = getattr(n3fit.hyper_optimization.rewards, hyper_loss)

        # Initialize the dictionaries which contain all fitting information
        self.input_list = []
        self.input_sizes = []
        self.training = {
            "output": [],
            "expdata": [],
            "losses": [],
            "ndata": 0,
            "model": None,
            "posdatasets": [],
            "posmultipliers": [],
            "posinitials": [],
            "integdatasets": [],
            "integmultipliers": [],
            "integinitials": [],
            "folds": [],
        }
        self.validation = {
            "output": [],
            "expdata": [],
            "losses": [],
            "ndata": 0,
            "model": None,
            "folds": [],
            "posdatasets": [],
        }
        self.experimental = {
            "output": [],
            "expdata": [],
            "losses": [],
            "ndata": 0,
            "model": None,
            "folds": [],
        }
        self.model_dicts = None

        self._fill_the_dictionaries()

        if self.validation["ndata"] == 0:
            # If there is no validation, the validation chi2 = training chi2
            self.no_validation = True
            self.validation["expdata"] = self.training["expdata"]
        else:
            # Consider the validation only if there is validation (of course)
            self.no_validation = False

        self.callbacks = []
        if debug:
            self.callbacks.append(callbacks.TimerCallback())

    def set_hyperopt(self, hyperopt_on, keys=None, status_ok="ok"):
        """ Set hyperopt options on and off (mostly suppresses some printing) """
        self.pass_status = status_ok
        if keys is None:
            keys = []
        self.hyperkeys = keys
        if hyperopt_on:
            self.print_summary = False
            self.mode_hyperopt = True
        else:
            self.print_summary = True
            self.mode_hyperopt = False

    ###########################################################################
    # # Internal functions                                                    #
    # Never to be called from the dark and cold outside world                 #
    ###########################################################################
    def _fill_the_dictionaries(self):
        """
        This function fills the following dictionaries
            -``training``: data for the fit
            -``validation``: data which for the stopping
            -``experimental``: 'true' data, only used for reporting purposes
        with fixed information.

        Fixed information: information which will not change between different runs of the code.
        This information does not depend on the parameters of the fit at any stage
        and so it will remain unchanged between different runs of the hyperoptimizer.

        The aforementioned information corresponds to:
            - ``expdata``: experimental data
            - ``name``: names of the experiment
            - ``ndata``: number of experimental points
        """
        for exp_dict in self.exp_info:
            self.training["expdata"].append(exp_dict["expdata"])
            self.validation["expdata"].append(exp_dict["expdata_vl"])
            self.experimental["expdata"].append(exp_dict["expdata_true"])

            self.training["folds"].append(exp_dict["folds"]["training"])
            self.validation["folds"].append(exp_dict["folds"]["validation"])
            self.experimental["folds"].append(exp_dict["folds"]["experimental"])

            nd_tr = exp_dict["ndata"]
            nd_vl = exp_dict["ndata_vl"]

            self.training["ndata"] += nd_tr
            self.validation["ndata"] += nd_vl
            self.experimental["ndata"] += nd_tr + nd_vl

            for dataset in exp_dict["datasets"]:
                self.all_datasets.append(dataset["name"])
        self.all_datasets = set(self.all_datasets)

        for pos_dict in self.pos_info:
            self.training["expdata"].append(pos_dict["expdata"])
            self.training["posdatasets"].append(pos_dict["name"])
            self.validation["expdata"].append(pos_dict["expdata"])
            self.validation["posdatasets"].append(pos_dict["name"])
        if self.integ_info is not None:
            for integ_dict in self.integ_info:
                self.training["expdata"].append(integ_dict["expdata"])
                self.training["integdatasets"].append(integ_dict["name"])

    def _model_generation(self, pdf_model, partition):
        """
        Fills the three dictionaries (``training``, ``validation``, ``experimental``)
        with the ``model`` entry


        Compiles the validation and experimental models with fakes optimizers and learning rate
        as they are never trained, but this is needed by some backends
        in order to run evaluate on them

        Before entering this function the dictionaries contain a list of inputs
        and a list of outputs, but they are not connected.
        This function connects inputs with outputs by injecting the PDF.
        At this point we have a PDF model that takes an input (1, None, 1)
        and outputs in return (1, none, 14)

        The injection of the PDF is done by concatenating all inputs and calling
        pdf_model on it.
        This in turn generates an output_layer that needs to be splitted for every experiment
        as we have a set of observable "functions" that each take (1, exp_xgrid_size, 14)
        and output (1, masked_ndata) where masked_ndata can be the training/validation
        or the experimental mask (in which cased masked_ndata == ndata)

        Parameters
        ----------
            pdf_model: MetaModel
                model producing pdf values

        Returns
        -------
            models: dict
                dict of MetaModels for training, validation and experimental
        """
        log.info("Generating the Model")

        # Compute the input array that will be given to the pdf
        input_arr = np.concatenate(self.input_list, axis=1)
        input_layer = operations.numpy_to_input(input_arr.T)

        # The input to the full model is expected to be the input to the PDF
        # by reutilizing `pdf_model.parse_input` we ensure any auxiliary input is also accunted fro
        full_model_input_dict = pdf_model._parse_input([input_layer], pass_content=False)

        # The output of the pdf on input_layer will be thus a concatenation
        # of the PDF values for all experiments
        full_pdf = pdf_model.apply_as_layer([input_layer])
        # The input layer is a concatenation of all experiments
        # we need now to split the output on a different array per experiment
        sp_ar = [self.input_sizes]
        sp_kw = {"axis": 1}
        splitting_layer = operations.as_layer(
            operations.split, op_args=sp_ar, op_kwargs=sp_kw, name="pdf_split"
        )
        splitted_pdf = splitting_layer(full_pdf)

        # If we are in a kfolding partition, select which datasets are out
        if partition:
            kfold_datasets = partition["datasets"]
            negate_k_datasets = [d for d in self.all_datasets if d not in kfold_datasets]
        else:
            kfold_datasets = None
            negate_k_datasets = None

        # Training and validation leave out the kofld dataset
        # experiment leaves out the negation
        output_tr = _pdf_injection(splitted_pdf, self.training["output"], kfold_datasets)
        training = MetaModel(full_model_input_dict, output_tr)

        output_vl = _pdf_injection(splitted_pdf, self.validation["output"], kfold_datasets)
        validation = MetaModel(full_model_input_dict, output_vl)

        output_ex = _pdf_injection(splitted_pdf, self.experimental["output"], negate_k_datasets)

        experimental = MetaModel(full_model_input_dict, output_ex)

        if self.print_summary:
            training.summary()

        models = {
            "training": training,
            "validation": validation,
            "experimental": experimental,
        }

        return models

    def _reset_observables(self):
        """
        Resets the 'output' and 'losses' entries of all 3 dictionaries:
                            (``training``, ``validation``, ``experimental``)
        as well as the input_list
        this is necessary as these can either depend on the parametrization of the NN
        or be obliterated when/if the backend state is reset
        """
        self.input_list = []
        self.input_sizes = []
        for key in ["output", "losses", "posmultipliers", "integmultipliers"]:
            self.training[key] = []
            self.validation[key] = []
            self.experimental[key] = []

    ############################################################################
    # # Parametizable functions                                                #
    #                                                                          #
    # The functions defined in this block accept a 'params' dictionary which   #
    # defines the fit and the behaviours of the Neural Networks                #
    #                                                                          #
    # These are all called by the function hyperparamizable below              #
    # i.e., the most important function is hyperparametrizable, which is a      #
    # wrapper around all of these                                              #
    ############################################################################
    def _generate_observables(
        self, all_pos_multiplier, all_pos_initial, all_integ_multiplier, all_integ_initial, epochs
    ):
        """
        This functions fills the 3 dictionaries (training, validation, experimental)
        with the output layers and the loss functions
        It also fill the list of input tensors (input_list)

        The arguments of this function are used to define the initial positivity of the
        positivity observables and the multiplier to be applied at each step.

        Parameters
        ----------
            all_pos_multiplier: float, None
                multiplier to be applied to the positivity each ``PUSH_POSITIVITY_EACH`` epochs
            all_pos_initial: float, None
                initial value for the positivity lambda
            epochs: int
                total number of epochs for the run
        """

        # First reset the dictionaries
        self._reset_observables()
        log.info("Generating layers")

        # Now we need to loop over all dictionaries (First exp_info, then pos_info and integ_info)
        for exp_dict in self.exp_info:
            if not self.mode_hyperopt:
                log.info("Generating layers for experiment %s", exp_dict["name"])

            exp_layer = model_gen.observable_generator(exp_dict)

            # Save the input(s) corresponding to this experiment
            self.input_list += exp_layer["inputs"]
            self.input_sizes.append(exp_layer["experiment_xsize"])

            # Now save the observable layer, the losses and the experimental data
            self.training["output"].append(exp_layer["output_tr"])
            self.validation["output"].append(exp_layer["output_vl"])
            self.experimental["output"].append(exp_layer["output"])

            self.training["losses"].append(exp_layer["loss_tr"])
            self.validation["losses"].append(exp_layer["loss_vl"])
            self.experimental["losses"].append(exp_layer["loss"])

        # Generate the positivity penalty
        for pos_dict in self.pos_info:
            if not self.mode_hyperopt:
                log.info("Generating positivity penalty for %s", pos_dict["name"])

            positivity_steps = int(epochs / PUSH_POSITIVITY_EACH)
            max_lambda = pos_dict["lambda"]

            pos_initial, pos_multiplier = _LM_initial_and_multiplier(
                all_pos_initial, all_pos_multiplier, max_lambda, positivity_steps
            )

            pos_layer = model_gen.observable_generator(pos_dict, positivity_initial=pos_initial)
            # The input list is still common
            self.input_list += pos_layer["inputs"]
            self.input_sizes.append(pos_layer["experiment_xsize"])

            # The positivity should be on both training and validation models
            self.training["output"].append(pos_layer["output_tr"])
            self.training["losses"].append(pos_layer["loss_tr"])
            self.validation["output"].append(pos_layer["output_tr"])
            self.validation["losses"].append(pos_layer["loss_tr"])

            self.training["posmultipliers"].append(pos_multiplier)
            self.training["posinitials"].append(pos_initial)

        # Finally generate the integrability penalty
        if self.integ_info is not None:
            for integ_dict in self.integ_info:
                if not self.mode_hyperopt:
                    log.info("Generating integrability penalty for %s", integ_dict["name"])

                integrability_steps = int(epochs / PUSH_INTEGRABILITY_EACH)
                max_lambda = integ_dict["lambda"]

                integ_initial, integ_multiplier = _LM_initial_and_multiplier(
                    all_integ_initial, all_integ_multiplier, max_lambda, integrability_steps
                )

                integ_layer = model_gen.observable_generator(
                    integ_dict, positivity_initial=integ_initial, integrability=True
                )
                # The input list is still common
                self.input_list += integ_layer["inputs"]
                self.input_sizes.append(integ_layer["experiment_xsize"])

                # The integrability all falls to the training
                self.training["output"].append(integ_layer["output_tr"])
                self.training["losses"].append(integ_layer["loss_tr"])
                self.training["integmultipliers"].append(integ_multiplier)
                self.training["integinitials"].append(integ_initial)

    def _generate_pdf(
        self,
        nodes_per_layer,
        activation_per_layer,
        initializer,
        layer_type,
        dropout,
        regularizer,
        regularizer_args,
        seed,
    ):
        """
        Defines the internal variable layer_pdf
        this layer takes any input (x) and returns the pdf value for that x

        if the sumrule is being imposed, it also updates input_list with the
        integrator_input tensor used to calculate the sumrule

        Parameters:
        -----------
            nodes_per_layer: list
                list of nodes each layer has
            activation_per_layer: list
                list of the activation function for each layer
            initializer: str
                initializer for the weights of the NN
            layer_type: str
                type of layer to be used
            dropout: float
                dropout to add at the end of the NN
            regularizer: str
                choice of regularizer to add to the dense layers of the NN
            regularizer_args: dict
                dictionary of arguments for the regularizer
            seed: int
                seed for the NN
        see model_gen.pdfNN_layer_generator for more information

        Returns
        -------
            pdf_model: MetaModel
                pdf model
        """
        log.info("Generating PDF model")

        # Set the parameters of the NN
        # Generate the NN layers
        pdf_model = model_gen.pdfNN_layer_generator(
            nodes=nodes_per_layer,
            activations=activation_per_layer,
            layer_type=layer_type,
            flav_info=self.flavinfo,
            fitbasis=self.fitbasis,
            seed=seed,
            initializer_name=initializer,
            dropout=dropout,
            regularizer=regularizer,
            regularizer_args=regularizer_args,
            impose_sumrule=self.impose_sumrule,
        )
        return pdf_model

    def _assign_data(self, models, fold_k=0):
        """Assign to each model the data to compare with as well as the
        number of data points in the model.
        In the most general case training and validation get assigned the replic'd data
        while experimental gets the actual data.
        In the kfolding case (i.e, partition != None), they all receive the same data
        but the folded data is set to 0 for training and validation
        """
        training = _assign_data_to_model(models["training"], self.training, fold_k)
        validation = _assign_data_to_model(models["validation"], self.validation, fold_k)
        experimental = _assign_data_to_model(models["experimental"], self.experimental, fold_k)

        ret = {
            "training": training,
            "validation": validation,
            "experimental": experimental,
        }
        return ret

    def _prepare_reporting(self, partition):
        """Parses the information received by the :py:class:`n3fit.ModelTrainer.ModelTrainer`
        to select the bits necessary for reporting the chi2.
        Receives the chi2 partition data to see whether any dataset is to be left out
        """
        reported_keys = ["name", "count_chi2", "positivity", "integrability", "ndata", "ndata_vl"]
        reporting_list = []
        for exp_dict in self.all_info:
            reporting_dict = {k: exp_dict.get(k) for k in reported_keys}
            if partition:
                # If we are in a partition we need to remove the number of datapoints
                # in order to avoid calculating the chi2 wrong
                for dataset in exp_dict["datasets"]:
                    if dataset in partition["datasets"]:
                        ndata = dataset["ndata"]
                        frac = dataset["frac"]
                        reporting_dict["ndata"] -= int(ndata * frac)
                        reporting_dict["ndata_vl"] = int(ndata * (1 - frac))
            reporting_list.append(reporting_dict)
        return reporting_list

    def _train_and_fit(self, training_model, stopping_object, epochs=100):
        """
        Trains the NN for the number of epochs given using
        stopping_object as the stopping criteria

        Every ``PUSH_POSITIVITY_EACH`` epochs the positivity will be multiplied by their
        respective positivity multipliers.
        In the same way, every ``PUSH_INTEGRABILITY_EACH`` epochs the integrability
        will be multiplied by their respective integrability multipliers
        """
        callback_st = callbacks.StoppingCallback(stopping_object)
        callback_pos = callbacks.LagrangeCallback(
            self.training["posdatasets"],
            self.training["posmultipliers"],
            update_freq=PUSH_POSITIVITY_EACH,
        )
        callback_integ = callbacks.LagrangeCallback(
            self.training["integdatasets"],
            self.training["integmultipliers"],
            update_freq=PUSH_INTEGRABILITY_EACH,
        )

        training_model.perform_fit(
            epochs=epochs,
            verbose=False,
            callbacks=self.callbacks + [callback_st, callback_pos, callback_integ],
        )

        if stopping_object.positivity_pass():
            return self.pass_status
        else:
            return self.failed_status

    def _hyperopt_override(self, params):
        """ Unrolls complicated hyperopt structures into very simple dictionaries"""
        # I love the smell of napalm in the morning
        for hyperkey in self.hyperkeys:
            item = params[hyperkey]
            if isinstance(item, dict):
                for key, value in item.items():
                    params[key] = value

    def enable_tensorboard(self, logdir, weight_freq = 0, profiling=False):
        """ Enables tensorboard callback for further runs of the fitting procedure 

        Parameters
        ----------
            logdir: Path
                path where to save the tensorboard logs
            weight_freq: int
                frequency (in epochs) at which to save weight histograms
            profiling: bool
                flag to enable the tensorboard profiler
        """
        callback_tb = callbacks.gen_tensorboard_callback(logdir, profiling=profiling, histogram_freq=weight_freq)
        self.callbacks.append(callback_tb)

    def evaluate(self, stopping_object):
        """Returns the training, validation and experimental chi2

        Parameters
        ----------
            stopping_object
                A Stopping intance which will have associated a validation model and the
                list of output layers that should contribute to the training chi2

        Returns
        -------
            train_chi2: chi2 of the trainining set
            val_chi2 : chi2 of the validation set
            exp_chi2: chi2 of the experimental data (without replica or tr/vl split)
        """
        if self.model_dicts is None:
            raise RuntimeError("Modeltrainer.evaluate was called before any training")
        # Needs to receive a `stopping_object` in order to select the part of the
        # training and the validation which are actually `chi2` and not part of the penalty
        training = self.model_dicts["training"]
        experimental = self.model_dicts["experimental"]
        train_chi2 = stopping_object.evaluate_training(training["model"])
        val_chi2, _ = stopping_object.validation.loss()
        exp_chi2 = experimental["model"].compute_losses()["loss"] / experimental["ndata"]
        return train_chi2, val_chi2, exp_chi2

    def hyperparametrizable(self, params):
        """
        Wrapper around all the functions defining the fit.

        After the ModelTrainer class has been instantiated,
        a call to this function (with a ``params`` dictionary) is necessary
        in order to generate the whole PDF model and perform a fit.

        This is a necessary step for hyperopt to work

        Parameters used only here:
            - ``epochs``: maximum number of iterations for the fit to run
            - ``stopping_patience``: patience of the stopper after finding a new minimum
        All other parameters are passed to the corresponding functions
        """

        # Reset the internal state of the backend every time this function is called
        print("")
        clear_backend_state()

        # Preprocess some hyperparameters
        epochs = int(params["epochs"])
        stopping_patience = params["stopping_patience"]
        stopping_epochs = int(epochs * stopping_patience)

        # When doing hyperopt some entries in the params dictionary
        # can bring with them overriding arguments
        if self.mode_hyperopt:
            log.info("Performing hyperparameter scan")
            for key in self.hyperkeys:
                log.info(" > > Testing %s = %s", key, params[key])
            self._hyperopt_override(params)

        # Fill the 3 dictionaries (training, validation, experimental) with the layers and losses
        # when k-folding, these are the same for all folds
        positivity_dict = params.get("positivity", {})
        integrability_dict = params.get("integrability", {})
        self._generate_observables(
            positivity_dict.get("multiplier"),
            positivity_dict.get("initial"),
            integrability_dict.get("multiplier"),
            integrability_dict.get("initial"),
            epochs,
        )
        threshold_pos = positivity_dict.get("threshold", 1e-6)
        threshold_chi2 = params.get("threshold_chi2", CHI2_THRESHOLD)

        # Initialize the chi2 dictionaries
        l_train = []
        l_valid = []
        l_exper = []
        l_hyper = []

        ### Training loop
        for k, partition in enumerate(self.kpartitions):
            # Each partition of the kfolding needs to have its own separate model
            seed = self.NNseed
            if k > 0:
                # Update the seed
                seed = np.random.randint(0, pow(2, 31))

            # Generate the pdf model
            pdf_model = self._generate_pdf(
                params["nodes_per_layer"],
                params["activation_per_layer"],
                params["initializer"],
                params["layer_type"],
                params["dropout"],
                params.get("regularizer", None),  # regularizer optional
                params.get("regularizer_args", None),
                seed,
            )

            # Model generation joins all the different observable layers
            # together with pdf model generated above
            models = self._model_generation(pdf_model, partition)

            # Only after model generation, apply possible weight file
            if self.model_file:
                log.info("Applying model file %s", self.model_file)
                pdf_model.load_weights(self.model_file)

            if k > 0:
                # Reset the positivity and integrability multipliers
                pos_and_int = self.training["posdatasets"] + self.training["integdatasets"]
                initial_values = self.training["posinitials"] + self.training["posinitials"]
                models["training"].reset_layer_weights_to(pos_and_int, initial_values)

            # Assign data to each model
            # model dicts is similar to model but includes information about
            # the target data and number of points
            model_dicts = self._assign_data(models, k)

            # Generate the list containing reporting info necessary for chi2
            reporting = self._prepare_reporting(partition)

            if self.no_validation:
                # Substitute the validation model with the training model
                model_dicts["validation"] = model_dicts["training"]
                validation_model = models["training"]
            else:
                validation_model = models["validation"]

            # Generate the stopping_object this object holds statistical information about the fit
            # it is used to perform stopping
            stopping_object = Stopping(
                validation_model,
                reporting,
                total_epochs=epochs,
                stopping_patience=stopping_epochs,
                save_weights_each=self.save_weights_each,
                threshold_positivity=threshold_pos,
                threshold_chi2=threshold_chi2
            )

            # Compile each of the models with the right parameters
            _model_compilation(model_dicts, params["optimizer"])

            passed = self._train_and_fit(
                models["training"],
                stopping_object,
                epochs=epochs,
            )

            # Compute validation and training loss
            training_loss = stopping_object.tr_loss
            validation_loss = stopping_object.vl_loss

            # Compute experimental loss
            exp_loss_raw = models["experimental"].compute_losses()["loss"]
            experimental_loss = exp_loss_raw / model_dicts["experimental"]["ndata"]

            if self.mode_hyperopt:
                hyper_loss = experimental_loss
                if passed != self.pass_status:
                    log.info("Hyperparameter combination fail to find a good fit, breaking")
                    # If the fit failed to fit, no need to add a penalty to the loss
                    break
                for penalty in self.hyper_penalties:
                    hyper_loss += penalty(pdf_model, stopping_object)
                l_hyper.append(hyper_loss)
                log.info("Fold %d finished, loss=%.1f, pass=%s", k+1, hyper_loss, passed)
                if hyper_loss > self.hyper_threshold:
                    log.info("Loss above threshold (%.1f > %.1f), breaking", hyper_loss, self.hyper_threshold)
                    # Apply a penalty proportional to the number of folds that have not been computed
                    pen_mul = len(self.kpartitions) - k
                    l_hyper = [i * pen_mul for i in l_hyper]
                    break

            # Save all losses
            l_train.append(training_loss)
            l_valid.append(validation_loss)
            l_exper.append(experimental_loss)

        dict_out = {
            "status": passed,
            "training_loss": np.average(l_train),
            "validation_loss": np.average(l_valid),
            "experimental_loss": np.average(l_exper),
        }

        if self.mode_hyperopt:
            dict_out["loss"] = self.hyper_loss(l_hyper)
            dict_out["kfold_meta"] = {
                "training_losses": l_train,
                "validation_losses": l_valid,
                "experimental_losses": l_exper,
                "hyper_losses": l_hyper,
            }
            # If we are using hyperopt we don't need to output any other information
            return dict_out

        dict_out["loss"] = experimental_loss

        # Add to the output dictionary things that are needed by performfit.py
        # to generate the output pdf, check the arc-length, gather stats, etc
        # some of them are already attributes of the class so they are redundant here
        # but I think it's good to present them explicitly
        dict_out["stopping_object"] = stopping_object
        dict_out["experimental"] = self.experimental
        dict_out["training"] = self.training
        dict_out["pdf_model"] = pdf_model

        # Only after the training has finished, we save all models for future reporting
        self.model_dicts = model_dicts

        return dict_out
