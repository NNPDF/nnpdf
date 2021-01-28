"""
    The ModelTrainer class is the true driver around the n3fit code

    This class is initialized with all information about the NN, inputs and outputs.
    The construction of the NN and the fitting is performed at the same time when the
    hyperparametrizable method of the function is called.

    This allows to use hyperscanning libraries, that need to change the parameters of the network
    between iterations while at the same time keeping the amount of redundant calls to a minimum
"""
import logging
from itertools import zip_longest
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


def _pdf_injection(pdf_layers, observables, masks):
    """
    Takes as input a list of output layers and returns a corresponding list
    where all output layers call the pdf layer at self.pdf_layer
    """
    return [f(x, mask=m) for f, x, m in zip_longest(observables, pdf_layers, masks)]


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
        kfold_parameters=None,
        max_cores=None,
        model_file=None,
        sum_rules=True,
        parallel_models=1
    ):
        """
        Parameters
        ----------
            exp_info: list
                list of dictionaries containing experiments
            pos_info: list
                list of dictionaries containing positivity sets
            integ_info: list
                list of dictionaries containing integrability sets
            flavinfo: list
                the object returned by fitting['basis']
            fitbasis: str
                the name of the basis being fitted
            nnseed: int
                the seed used to initialise the Neural Network, will be passed to model_gen
            pass_status: str
                flag to signal a good run
            failed_status: str
                flag to signal a bad run
            debug: bool
                flag to activate some debug options
            kfold_parameters: dict
                parameters defining the kfolding method
            max_cores: int
                maximum number of cores the fitting can use to run
            model_file: str
                whether to save the models
            sum_rules: str
		whether sum rules should be enabled (All, MSR, VSR, False)
            parallel_models: int
                number of models to fit in parallel
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
        self.all_datasets = []
        self.parallel_models = parallel_models

        # Initialise internal variables which define behaviour
        if debug:
            self.max_cores = 1
        else:
            self.max_cores = max_cores
        self.model_file = model_file
        self.print_summary = True
        self.mode_hyperopt = False
        self.impose_sumrule = sum_rules
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
            "ndata": 0,
            "model": None,
            "folds": [],
            "posdatasets": [],
        }
        self.experimental = {
            "output": [],
            "expdata": [],
            "ndata": 0,
            "model": None,
            "folds": [],
        }
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

    def _model_generation(self, pdf_models, partition, partition_idx):
        """
        Fills the three dictionaries (``training``, ``validation``, ``experimental``)
        with the ``model`` entry.

        Compiles the validation and experimental models with fakes optimizers and learning rate
        as they are never trained, but this is needed by some backends
        in order to run evaluate on them.

        Before entering this function the dictionaries contain a list of inputs
        and a list of outputs, but they are not connected.
        This function connects inputs with outputs by injecting the PDF.
        At this point we have a PDF model that takes an input (1, None, 1)
        and outputs in return (1, none, 14).

        The injection of the PDF is done by concatenating all inputs and calling
        pdf_model on it.
        This in turn generates an output_layer that needs to be splitted for every experiment
        as we have a set of observable "functions" that each take (1, exp_xgrid_size, 14)
        and output (1, masked_ndata) where masked_ndata can be the training/validation
        or the experimental mask (in which cased masked_ndata == ndata).

        Several models can be fitted at once by passing a list of models with a shared input
        this function will give the same input to every model and will concatenate the output at the end
        so that the final output of the model is (1, None, 14, n) (with n=number of parallel models)


        Parameters
        ----------
            pdf_models: list(n3fit.backend.MetaModel)
                a list of models that produce PDF values

        Returns
        -------
            models: dict
                dict of MetaModels for training, validation and experimental
        """
        log.info("Generating the Model")

        # Compute the input array that will be given to the pdf
        input_arr = np.concatenate(self.input_list, axis=1)
        input_layer = operations.numpy_to_input(input_arr.T)

        # The trainable part of the model is a concatenation of all PDF models
        # where each model corresponds to a different replica
        all_replicas_pdf = []

        for pdf_model in pdf_models:
            # The input to the full model is expected to be the input to the PDF
            # by reutilizing `pdf_model.parse_input` we ensure any auxiliary input is also accunted fro
            full_model_input_dict = pdf_model._parse_input([input_layer], pass_content=False)

            # The output of the pdf on input_layer will be thus a concatenation
            # of the PDF values for all experiments
            full_pdf = pdf_model.apply_as_layer([input_layer])

            all_replicas_pdf.append(full_pdf)

        full_pdf_per_replica = operations.stack(all_replicas_pdf, axis=-1)

        # The input layer is a concatenation of all experiments
        # we need now to split the output on a different array per experiment
        sp_ar = [self.input_sizes]
        sp_kw = {"axis": 1}
        splitting_layer = operations.as_layer(
            operations.split, op_args=sp_ar, op_kwargs=sp_kw, name="pdf_split"
        )
        splitted_pdf = splitting_layer(full_pdf_per_replica)

        # If we are in a kfolding partition, select which datasets are out
        if partition:
            training_mask = [i[partition_idx] for i in self.training["folds"]]
            validation_mask = [i[partition_idx] for i in self.validation["folds"]]
            experimental_mask = [i[partition_idx] for i in self.experimental["folds"]]
        else:
            training_mask = validation_mask = experimental_mask = [None]

        # Training and validation leave out the kofld dataset
        # experiment leaves out the negation
        output_tr = _pdf_injection(splitted_pdf, self.training["output"], training_mask)
        training = MetaModel(full_model_input_dict, output_tr)

        # For validation we don't have integrability
        n_val = len(self.validation["output"])
        output_vl = _pdf_injection(splitted_pdf[:n_val], self.validation["output"], validation_mask)
        validation = MetaModel(full_model_input_dict, output_vl)

        # And for experimental we don't have positivity
        # TODO: risky to rely that much in the order
        n_exps = len(self.experimental["output"])
        output_ex = _pdf_injection(
            splitted_pdf[:n_exps], self.experimental["output"], experimental_mask
        )

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
        for key in ["output", "posmultipliers", "integmultipliers"]:
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
            self.validation["output"].append(pos_layer["output_tr"])

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
        log.info("Generating PDF models")

        # Set the parameters of the NN
        # Generate the NN layers
        pdf_models = model_gen.pdfNN_layer_generator(
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
            parallel_models=self.parallel_models,
        )
        return pdf_models

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

        # TODO: in order to use multireplica in hyperopt is is necessary to define what "passing" means
        # for now consider the run as good if any replica passed
        if any([bool(i) for i in stopping_object.e_best_chi2]):
            return self.pass_status
        return self.failed_status

    def _hyperopt_override(self, params):
        """ Unrolls complicated hyperopt structures into very simple dictionaries"""
        # I love the smell of napalm in the morning
        for hyperkey in self.hyperkeys:
            item = params[hyperkey]
            if isinstance(item, dict):
                for key, value in item.items():
                    params[key] = value

    def enable_tensorboard(self, logdir, weight_freq=0, profiling=False):
        """Enables tensorboard callback for further runs of the fitting procedure

        Parameters
        ----------
            logdir: Path
                path where to save the tensorboard logs
            weight_freq: int
                frequency (in epochs) at which to save weight histograms
            profiling: bool
                flag to enable the tensorboard profiler
        """
        callback_tb = callbacks.gen_tensorboard_callback(
            logdir, profiling=profiling, histogram_freq=weight_freq
        )
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
        if self.training["model"] is None:
            raise RuntimeError("Modeltrainer.evaluate was called before any training")
        # Needs to receive a `stopping_object` in order to select the part of the
        # training and the validation which are actually `chi2` and not part of the penalty
        train_chi2 = stopping_object.evaluate_training(self.training["model"])
        val_chi2 = stopping_object.vl_chi2
        exp_chi2 = self.experimental["model"].compute_losses()["loss"] / self.experimental["ndata"]
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
            pdf_models = self._generate_pdf(
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
            models = self._model_generation(pdf_models, partition, k)

            # Only after model generation, apply possible weight file
            # TODO: not sure whether it is a good idea that all of them start at the same point
            if self.model_file:
                log.info("Applying model file %s", self.model_file)
                for pdf_model in pdf_models:
                    pdf_model.load_weights(self.model_file)

            if k > 0:
                # Reset the positivity and integrability multipliers
                pos_and_int = self.training["posdatasets"] + self.training["integdatasets"]
                initial_values = self.training["posinitials"] + self.training["posinitials"]
                models["training"].reset_layer_weights_to(pos_and_int, initial_values)

            # Generate the list containing reporting info necessary for chi2
            reporting = self._prepare_reporting(partition)

            if self.no_validation:
                # Substitute the validation model with the training model
                models["validation"] = models["training"]
                validation_model = models["training"]
            else:
                validation_model = models["validation"]

            # Generate the stopping_object this object holds statistical information about the fit
            # it is used to perform stopping
            stopping_object = Stopping(
                validation_model,
                reporting,
                pdf_models,
                total_epochs=epochs,
                stopping_patience=stopping_epochs,
                threshold_positivity=threshold_pos,
                threshold_chi2=threshold_chi2,
            )

            # Compile each of the models with the right parameters
            for model in models.values():
                model.compile(**params["optimizer"])

            passed = self._train_and_fit(
                models["training"],
                stopping_object,
                epochs=epochs,
            )

            if self.mode_hyperopt:
                # TODO: currently only working for one single replica
                # If doing a hyperparameter scan we need to keep track at this point of the loss function
                validation_loss = stopping_object.vl_chi2

                # Compute experimental loss
                exp_loss_raw = np.take(models["experimental"].compute_losses()["loss"], -1)
                # And divide by the number of active points in this fold
                # it would be nice to have a ndata_per_fold variable coming in the vp object...
                ndata = np.sum([np.count_nonzero(i[k]) for i in self.experimental["folds"]])
                experimental_loss = exp_loss_raw / ndata

                hyper_loss = experimental_loss
                if passed != self.pass_status:
                    log.info("Hyperparameter combination fail to find a good fit, breaking")
                    # If the fit failed to fit, no need to add a penalty to the loss
                    break
                for penalty in self.hyper_penalties:
                    hyper_loss += penalty(pdf_model, stopping_object)
                l_hyper.append(hyper_loss)
                log.info("Fold %d finished, loss=%.1f, pass=%s", k + 1, hyper_loss, passed)
                if hyper_loss > self.hyper_threshold:
                    log.info(
                        "Loss above threshold (%.1f > %.1f), breaking",
                        hyper_loss,
                        self.hyper_threshold,
                    )
                    # Apply a penalty proportional to the number of folds that have not been computed
                    pen_mul = len(self.kpartitions) - k
                    l_hyper = [i * pen_mul for i in l_hyper]
                    break

                # Save all losses
                l_valid.append(validation_loss)
                l_exper.append(experimental_loss)

        if self.mode_hyperopt:
            # Hyperopt needs a dictionary with information about the losses
            # it is possible to store arbitrary information in the trial file by adding it to this dictionary
            dict_out = {
                "status": passed,
                "loss": self.hyper_loss(l_hyper),
                "validation_loss": np.average(l_valid),
                "experimental_loss": np.average(l_exper),
                "kfold_meta": {
                    "validation_losses": l_valid,
                    "experimental_losses": l_exper,
                    "hyper_losses": l_hyper,
                },
            }
            return dict_out

        # Keep a reference to the models after training for future reporting
        self.training["model"] = models["training"]
        self.experimental["model"] = models["experimental"]
        self.validation["model"] = models["validation"]

        # In a normal run, the only information we need to output is the stopping object
        # (which contains metadata about the stopping)
        # and the pdf models (which are used to generate the PDF grids and compute arclengths)
        dict_out = {"status": passed, "stopping_object": stopping_object, "pdf_models": pdf_models}
        return dict_out
