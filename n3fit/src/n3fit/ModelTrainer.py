"""
    The ModelTrainer class is the true driver around the n3fit code

    This class is initialized with all information about the NN, inputs and outputs.
    The construction of the NN and the fitting is performed at the same time when the
    hyperparametizable method of the function is called.

    This allows to use hyperscanning libraries, that need to change the parameters of the network
    between iterations while at the same time keeping the amount of redundant calls to a minimum
"""
import logging
import numpy as np
import n3fit.model_gen as model_gen
from n3fit.backends import MetaModel, clear_backend_state, operations
from n3fit.stopping import Stopping

log = logging.getLogger(__name__)

# Threshold defaults
# Any partition with a chi2 over the threshold will discard its hyperparameters
HYPER_THRESHOLD = 5.0
# The stopping will not consider any run where the validation is not under this threshold
THRESHOLD_CHI2 = 10.0


def _fold_data(all_data, folds, k_idx, negate_fold=False):
    """
    Utility method to fold the data.
    If the folds are emtpy returns the data unchganed.
    If the number of folds and data are different, it is assumed that anything beyond
    the fold passes through.
    If negate_fold is active the fold is negated before applying it.

    Parameters
    ----------
        all_data: list of np.array
            original data
        folds: list of (list of) np.array
            boolean array to select the data that goes through
        k_idx: int
            index of the fold
        negate_fold: bool
            Flag to decide whether the fold such be negated

    Returns
    -------
        folded_data: list of np.array
            subset of the original data
    """
    if any(folds):
        # Check how many datasets we can fold
        folded_data = []
        nfolds = len(folds)
        for original_data, fold in zip(all_data, folds):
            if negate_fold:
                kfold = ~fold[k_idx]
            else:
                kfold = fold[k_idx]
            folded_data.append(original_data * kfold)
        # If there are dataset from the original set, add them all
        folded_data += all_data[nfolds:]
        return folded_data
    else:
        return all_data


def _compile_one_model(model_dict, kidx=None, negate_fold=False, **params):
    """
    Compiles one model dictionary. The model dictionary must include a backend-dependent
    model (`model`), a list of losses (`losses`), data to be compared with (`data`) and,
    if applies, a "fold".

    Parameters
    ----------
        model_dict: dict
            A ditionary defining the model
        kidx: int
            k-index of the fold
        negate_fold: bool
            Flag to pass to `fold_data` to negate the fold
        **params: **dict
            Parameters to be passes to the compile method of the model
    """
    model = model_dict["model"]
    losses = model_dict["losses"]
    data = model_dict["expdata"]
    fold = model_dict["folds"]
    folded_data = _fold_data(data, fold, kidx, negate_fold=negate_fold)
    model.compile(loss=losses, target_output=folded_data, **params)


def _count_data_fold(folds, k_idx):
    """
    Count the number of active datapoints in this fold

    Parameters
    ----------
        folds: list
            list of folds
        k_idx: int
            index of the fold to count
    """
    n = 0
    for fold in folds:
        n += np.count_nonzero(~fold[k_idx])
    return n


def _pdf_injection(pdf_layers, observables):
    """
    Takes as input a list of output layers and returns a corresponding list
    where all output layers call the pdf layer at self.pdf_layer
    """
    return [f(x) for f, x in zip(observables, pdf_layers)]


class ModelTrainer:
    """
        ModelTrainer Class:

        Wrapper around the fitting code and the generation of the Neural Network

        When the "hyperparametizable"* function is called with a dictionary of parameters,
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
        flavinfo,
        nnseed,
        pass_status="ok",
        failed_status="fail",
        debug=False,
        save_weights_each=False,
        kfold_parameters=None,
        max_cores=None,
    ):
        """
        # Arguments:
            - `exp_info`: list of dictionaries containing experiments
            - `pos_info`: list of dictionaries containing positivity sets
            - `flavinfo`: the object returned by fitting['basis']
            - `nnseed`: the seed used to initialise the Neural Network, will be passed to model_gen
            - `pass_status`: flag to signal a good run
            - `failed_status`: flag to signal a bad run
            - `pass_status`: flag to signal a good run
            - `failed_status`: flag to signal a bad run
            - `debug`: flag to activate some debug options
            - `save_weights_each`: if set, save the state of the fit
                                    every `save_weights_each` epochs
        """

        # Save all input information
        self.exp_info = exp_info
        self.pos_info = pos_info
        self.all_info = exp_info + pos_info
        self.flavinfo = flavinfo
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
        self.print_summary = True
        self.mode_hyperopt = False
        self.model_file = None
        self.impose_sumrule = True
        self.hyperkeys = None
        if kfold_parameters is None:
            self.kpartitions = [None]
            self.hyper_threshold = None
        else:
            self.kpartitions = kfold_parameters["partitions"]
            self.hyper_threshold = kfold_parameters.get("threshold", HYPER_THRESHOLD)

        # Initialize the pdf model
        self.pdf_model = None
        self.integrator_input = None

        # Initialize the dictionaries which contain all fitting information
        self.input_list = []
        self.input_sizes = []
        self.training = {
            "output": [],
            "expdata": [],
            "losses": [],
            "original_ndata": 0,
            "ndata": 0,
            "model": None,
            "posdatasets": [],
            "folds": [],
        }
        self.validation = {
            "output": [],
            "expdata": [],
            "losses": [],
            "original_ndata": 0,
            "ndata": 0,
            "model": None,
            "folds": [],
        }
        self.experimental = {
            "output": [],
            "expdata": [],
            "losses": [],
            "original_ndata": 0,
            "ndata": 0,
            "model": None,
            "folds": [],
        }
        self.list_of_models_dicts = [self.training, self.experimental]

        self._fill_the_dictionaries()

        if self.validation["ndata"] == 0:
            # If there is no validation, the validation chi2 = training chi2
            self.no_validation = True
            self.validation["expdata"] = self.training["expdata"]
        else:
            # Consider the validation only if there is validation (of course)
            self.no_validation = False
            self.list_of_models_dicts.append(self.validation)

        # Save the original number of datapoints in case we need to overwrite it later
        for model_dict in self.list_of_models_dicts:
            model_dict["original_ndata"] = model_dict["ndata"]

    @property
    def model_file(self):
        """ If a model_file is set the training model will try to get the weights form here """
        return self._model_file

    @model_file.setter
    def model_file(self, model_file):
        self._model_file = model_file

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
            -`training`: data for the fit
            -`validation`: data which for the stopping
            -`experimental`: 'true' data, only used for reporting purposes
        with fixed information.

        Fixed information: information which will not change between different runs of the code.
        This information does not depend on the parameters of the fit at any stage
        and so it will remain unchanged between different runs of the hyperoptimizer.

        The aforementioned information corresponds to:
            - `expdata`: experimental data
            - `name`: names of the experiment
            - `ndata`: number of experimental points
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

    def _model_generation(self):
        """
        Fills the three dictionaries (`training`, `validation`, `experimental`)
        with the `model` entry


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
        """
        log.info("Generating the Model")

        # Compute the input array that will be given to the pdf
        input_arr = np.concatenate(self.input_list, axis=1)
        input_layer = operations.numpy_to_input(input_arr.T)
        if self.impose_sumrule:
            full_model_input = [self.integrator_input, input_layer]
        else:
            full_model_input = [input_layer]
        # The output of the pdf on input_layer will be thus a concatenation
        # of the PDF values for all experiments
        full_pdf = self.pdf_model.apply_as_layer([input_layer])
        # The input layer is a concatenation of all experiments
        # we need now to split the output on a different array per experiment
        sp_ar = [self.input_sizes]
        sp_kw = {"axis": 1}
        splitting_layer = operations.as_layer(
            operations.split, op_args=sp_ar, op_kwargs=sp_kw, name="pdf_split"
        )
        splitted_pdf = splitting_layer(full_pdf)

        for model_dict in self.list_of_models_dicts:
            output = _pdf_injection(splitted_pdf, model_dict["output"])
            model_dict["model"] = MetaModel(full_model_input, output)

        if self.model_file:
            # If a model file is given, load the weights from there
            # note: even though the load corresponds to the training model only,
            #       the layer_pdf is shared  and so it should affect all models
            self.training["model"].load_weights(self.model_file)

        if self.print_summary:
            self.training["model"].summary()

    def _reset_observables(self):
        """
        Resets the 'output' and 'losses' entries of all 3 dictionaries:
                            (`training`, `validation`, `experimental`)
        as well as the input_list
        this is necessary as these can either depend on the parametrization of the NN
        or be obliterated when/if the backend state is reset
        """
        self.input_list = []
        self.input_sizes = []
        for key in ["output", "losses"]:
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
    # i.e., the most important function is hyperparametizable, which is a      #
    # wrapper around all of these                                              #
    ############################################################################
    def _generate_observables(self, pos_multiplier, pos_initial):
        """
        This functions fills the 3 dictionaries (training, validation, experimental)
        with the output layers and the loss functions
        It also fill the list of input tensors (input_list)

        Parameters accepted:
            - `pos_multiplier`: the multiplier to be applied to the positivity each 100 epochs
            - `pos_initial`: the initial value for the positivity
        """

        # First reset the dictionaries
        self._reset_observables()
        log.info("Generating layers")

        # Now we need to loop over all dictionaries (First exp_info, then pos_info)
        for exp_dict in self.exp_info:
            if not self.mode_hyperopt:
                log.info("Generating layers for experiment %s", exp_dict["name"])

            exp_layer = model_gen.observable_generator(
                exp_dict, kfolding=self.mode_hyperopt
            )

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

        # Finally generate the positivity penalty
        for pos_dict in self.pos_info:
            if not self.mode_hyperopt:
                log.info("Generating positivity penalty for %s", pos_dict["name"])
            pos_layer = model_gen.observable_generator(
                pos_dict,
                positivity_initial=pos_initial,
                positivity_multiplier=pos_multiplier,
                kfolding=self.mode_hyperopt,
            )
            # The input list is still common
            self.input_list += pos_layer["inputs"]
            self.input_sizes.append(pos_layer["experiment_xsize"])

            # The positivity all falls to the training
            self.training["output"].append(pos_layer["output_tr"])
            self.training["losses"].append(pos_layer["loss_tr"])
        # Save the positivity multiplier into the training dictionary
        # as it will be used during training
        self.training["pos_multiplier"] = pos_multiplier

    def _generate_pdf(
        self,
        nodes_per_layer,
        activation_per_layer,
        initializer,
        layer_type,
        dropout,
        regularizer,
        regularizer_args,
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
        see model_gen.pdfNN_layer_generator for more information

        Returns
        -------
            layers: list
                list of layers
            integrator_input: tensor
                xgrid used for integration within the network

        note: both are being used at the moment for reporting purposes at the end of the fit
        """
        log.info("Generating PDF layer")

        # Set the parameters of the NN

        # Generate the NN layers
        pdf_model, integrator_input = model_gen.pdfNN_layer_generator(
            nodes=nodes_per_layer,
            activations=activation_per_layer,
            layer_type=layer_type,
            flav_info=self.flavinfo,
            seed=self.NNseed,
            initializer_name=initializer,
            dropout=dropout,
            regularizer=regularizer,
            regularizer_args=regularizer_args,
            impose_sumrule=self.impose_sumrule,
        )
        self.integrator_input = integrator_input
        self.pdf_model = pdf_model

    def _toggle_fold(self, datasets, kidx=0, off=False, recompile=False):
        """ Toggle the input dataset on (off) and turn all other datasets off (on)

        Parameters
        ----------
            datasets: list
                list of datasets defining the fold
            kidx: int
                index of the fold
            off: bool
                Choose whether to the input datasets should be on or off
            recompile: bool
                Choose whether the model should be recompiled
        """  # TODO datasets and kidx is redundant but...
        # TODO fold the number of datapoints in the dictionaries
        all_other_datasets = self.all_datasets - set(datasets)

        if off:
            val = 0.0
        else:
            val = 1.0

        # Use the "experimental" model to access all layers
        self.experimental["model"].set_masks_to(datasets, val=val)
        self.experimental["model"].set_masks_to(all_other_datasets, val=1.0 - val)

        # Now run through all dicts and count the actual total numer of points
        for model_dict in self.list_of_models_dicts:
            fold_ndata = _count_data_fold(model_dict["folds"], kidx)
            if off:
                new_ndata = model_dict["original_ndata"] - fold_ndata
            else:
                new_ndata = fold_ndata
            model_dict["ndata"] = new_ndata

        if recompile:
            _compile_one_model(self.experimental, kidx=kidx, negate_fold=not off)

    def _model_compilation(self, learning_rate, optimizer, kidx=None):
        """
        Wrapper around `_compile_one_model` to pass the right parameters
        and index of the k-folding.

        Currently the accepted for compilation are `learning_rate` and `optimizer`

        Parameters
        ----------
            learning_rate: float
                value of the learning rate
            optimizer: str
                name of the optimizer to be used
                optimizers accepted are backend-dependent
            kidx: int
                k-index of the fold
        """

        # Compile all different models
        for model_dict in self.list_of_models_dicts:
            param_dict = {"learning_rate": learning_rate, "optimizer_name": optimizer}
            _compile_one_model(model_dict, kidx=kidx, **param_dict)

    def _train_and_fit(self, stopping_object, epochs):
        """
        Trains the NN for the number of epochs given using
        stopping_object as the stopping criteria

        Every 100 epochs the positivitiy will be updated with
        self.training['pos_multiplier']
        """
        training_model = self.training["model"]
        pos_multiplier = self.training["pos_multiplier"]
        # Train the model for the number of epochs given
        for epoch in range(epochs):
            out = training_model.perform_fit(verbose=False)
            print_stats = False

            if (epoch + 1) % 100 == 0:
                print_stats = True
                training_model.multiply_weights(
                    self.training["posdatasets"], pos_multiplier
                )

            passes = stopping_object.monitor_chi2(out, epoch, print_stats=print_stats)

            if stopping_object.stop_here():
                break

        # Report a "good" training only if there was no NaNs
        # and if there was at least a point which passed positivity
        if passes and stopping_object.positivity:
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

    def evaluate(self, stopping_object):
        """ Returns the training, validation and experimental chi2

        Parameters
        ----------
            `stopping_object`
                A Stopping intance which will have associated a validation model and the
                list of output layers that should contribute to the training chi2

        Returns
        -------
            `train_chi2`: chi2 of the trainining set
            `val_chi2` : chi2 of the validation set
            `exp_chi2`: chi2 of the experimental data (without replica or tr/vl split)
        """
        # Needs to receive a `stopping_object` in order to select the part of the
        # training and the validation which are actually `chi2` and not part of the penalty
        train_chi2 = stopping_object.evaluate_training(self.training["model"])
        val_chi2, _ = stopping_object.validation.loss()
        exp_chi2 = (
            self.experimental["model"].compute_losses()["loss"]
            / self.experimental["ndata"]
        )
        return train_chi2, val_chi2, exp_chi2

    def hyperparametizable(self, params):
        """
        Wrapper around all the functions defining the fit.

        After the ModelTrainer class has been instantiated,
        a call to this function (with a `params` dictionary) is necessary
        in order to generate the whole PDF model and perform a fit.

        This is a necessary step for hyperopt to work

        Parameters used only here:
            - `epochs`: maximum number of iterations for the fit to run
            - `stopping_patience`: patience of the stopper after finding a new minimum
        All other parameters are passed to the corresponding functions
        """

        # Reset the internal state of the backend
        print("")
        if not self.debug or self.mode_hyperopt:
            clear_backend_state(max_cores=self.max_cores)

        # When doing hyperopt some entries in the params dictionary
        # can bring with them overriding arguments
        if self.mode_hyperopt:
            log.info("Performing hyperparameter scan")
            for key in self.hyperkeys:
                log.info(" > > Testing %s = %s", key, params[key])
            self._hyperopt_override(params)
            hyper_losses = []

        # Fill the 3 dictionaries (training, validation, experimental) with the layers and losses
        self._generate_observables(params["pos_multiplier"], params["pos_initial"])

        # Generate the pdf model
        self._generate_pdf(
            params["nodes_per_layer"],
            params["activation_per_layer"],
            params["initializer"],
            params["layer_type"],
            params["dropout"],
            params.get("regularizer", None),  # regularizer optional
            params.get("regularizer_args", None),
        )

        # Model generation
        self._model_generation()

        # Generate the stopping_object
        # this object holds statistical information about the fit
        # it can be used to perform stopping
        epochs = int(params["epochs"])
        stopping_patience = params["stopping_patience"]
        stopping_epochs = epochs * stopping_patience

        if self.no_validation:
            validation_model = self.training["model"]
        else:
            validation_model = self.validation["model"]

        # Initialize the chi2 dictionaries
        l_train = []
        l_valid = []
        l_exper = []

        ### Training loop
        for k, partition in enumerate(self.kpartitions):

            if self.mode_hyperopt:
                # Disable the partition for training
                datasets = partition["datasets"]
                self._toggle_fold(datasets, kidx=k, off=True)

            # Generate the stopping object
            stopping_object = Stopping(
                validation_model,
                self.all_info,
                total_epochs=epochs,
                stopping_patience=stopping_epochs,
                threshold_chi2=params.get("threshold_chi2", THRESHOLD_CHI2),
                save_weights_each=self.save_weights_each,
            )

            # Compile the training['model'] with the given parameters
            self._model_compilation(
                params["learning_rate"], params["optimizer"], kidx=k
            )

            passed = self._train_and_fit(stopping_object, epochs)

            # Compute validation and training loss
            training_loss = stopping_object.tr_loss
            validation_loss = stopping_object.vl_loss

            # Compute experimental loss
            exp_loss_raw = self.experimental["model"].compute_losses()["loss"]
            experimental_loss = exp_loss_raw / self.experimental["ndata"]

            l_train.append(training_loss)
            l_valid.append(validation_loss)
            l_exper.append(experimental_loss)

            if self.mode_hyperopt:
                # Toggle the fold
                self._toggle_fold(datasets, kidx=k, off=False, recompile=True)
                # Compute the hyperopt loss
                hyper_loss = (
                    self.experimental["model"].compute_losses()["loss"]
                    / self.experimental["ndata"]
                )
                hyper_losses.append(hyper_loss)
                # Check whether this run is any good, if not, get out
                log.info(f"fold: {k}")
                log.info(f"Experimental loss: {experimental_loss}")
                log.info(f"Hyper loss: {hyper_loss}")
                # Check whether the losses are above the set thresholds
                all_losses = [training_loss, validation_loss, experimental_loss]
                if any([il > self.hyper_threshold for il in all_losses]):
                    # Add a penalty for leaving early
                    hyper_losses = [ih + self.hyper_threshold for ih in hyper_losses]
                    log.info("Bad threshold: break")
                    break
                self.training["model"].reinitialize()

        dict_out = {
            "status": passed,
            "training_loss": np.average(l_train),
            "validation_loss": np.average(l_valid),
            "experimental_loss": np.average(l_exper),
        }

        if self.mode_hyperopt:
            ave = np.average(hyper_losses)
            std = np.var(hyper_losses)
            dict_out["loss"] = ave
            dict_out["kfold_meta"] = {
                "training_losses": l_train,
                "validation_losses": l_valid,
                "experimental_losses": l_exper,
                "hyper_losses": hyper_losses,
                "hyper_avg": ave,
                "hyper_std": std,
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
        dict_out["pdf_model"] = self.pdf_model

        return dict_out
