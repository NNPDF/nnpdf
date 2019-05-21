"""
    ModelTrainer Class:

    This class provides a wrapper around the fitting code and the generation of the Neural Network
    When the "hyperparametizable"* function is called with a dictionary of parameters,
    it generates a NN and subsequentially performs a fit.

    The motivation behind this class is minimising the amount of redundant calls of each hyperopt
    run, in particular this allows to completely reset the NN at the beginning of each iteration
    reusing some of the previous work.

    *called in this way because it accept a dictionary of parameters which defines the Neural Network
"""
import model_gen
import constrains
import Statistics
from backends import MetaModel, clear_backend_state

# If an experiment matches this name, it will not be included int he training
TESTNAME = "TEST"
# Weights for the hyperopt cost function
test_multiplier = 0.5
validation_multiplier = 0.5


class ModelTrainer:
    def __init__(
        self, exp_info, pos_info, flavinfo, nnseed, pass_status="ok", failed_status="fail", log=None, debug=False
    ):
        """
        - exp_info: list of dictionaries containing experiments
        - pos_info: list of dictionaries containing positivity sets
        - flavinfo: the object returned by fitting['basis']
        - nnseed: the seed used to initialise the Neural Network, will be passed to model_gen
        - log: if given, printing will be done to log.info
        - debug: flag to activate some debug options
        """

        # Save all input information
        self.exp_info = exp_info
        self.pos_info = pos_info
        self.flavinfo = flavinfo
        self.NNseed = nnseed
        self.pass_status = pass_status
        self.failed_status = failed_status
        self.debug = debug

        # Initialise internal variables which define behaviour
        self.print_summary = True
        self.mode_hyperopt = False
        self.model_file = None
        self.impose_sumrule = True

        # Set the log_info function
        if log:
            self.log_info = log.info
        else:
            self.log_info = print
        self.log = log

        # Initialize the pdf layer
        self.layer_pdf = None

        # Initialize the dictionaries which contain all fitting information
        self.input_list = []
        self.ndata_dict = {}
        self.training = {"output": [], "expdata": [], "losses": [], "ndata": 0, "model": None, "posdatasets": []}
        self.validation = {"output": [], "expdata": [], "losses": [], "ndata": 0, "model": None}
        self.experimental = {"output": [], "expdata": [], "losses": [], "ndata": 0, "model": None}

        self.test_dict = None
        self._fill_the_dictionaries()

        if self.validation["ndata"] == 0:
            self.no_validation = True
            new_dict = {}
            for key, item in self.ndata_dict.items():
                if key.endswith("_vl"):
                    continue
                new_dict[key + "_vl"] = item
            self.ndata_dict.update(new_dict)
            self.validation["expdata"] = self.training["expdata"]
        else:
            self.no_validation = False

    def set_model_file(self, model_file):
        """ If a model_file is set the training model will try to get the weights form here """
        self.model_file = model_file

    def set_hyperopt(self, on, keys=None):
        """ Set hyperopt options on and off (mostly suppresses some printing) """
        if keys is None:
            keys = []
        self.hyperkeys = keys
        if on:
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
        This function fills the following dictionaries with fixed information,
            - training: data for the fit
            - validation: data which for the stopping
            - experimental: 'true' data, only used for reporting purposes
            - ndata_dict: a dictionary containing 'name of experiment' : Number Of Points

        the fixed information (i.e., non dependant on the parameters of the fit) is:
            - experimental data (expdata)
            - names of the experiment (name)
            - number of experimental points (ndata)
        """
        for exp_dict in self.exp_info:
            if exp_dict["name"] == TESTNAME:
                self.test_dict = {
                    "output": [],
                    "expdata": exp_dict["expdata_true"],
                    "losses": [],
                    "ndata": exp_dict["ndata"] + exp_dict["ndata_vl"],
                    "model": None,
                }
                continue
            self.training["expdata"].append(exp_dict["expdata"])
            self.validation["expdata"].append(exp_dict["expdata_vl"])
            self.experimental["expdata"].append(exp_dict["expdata_true"])

            nd_tr = exp_dict["ndata"]
            nd_vl = exp_dict["ndata_vl"]

            self.ndata_dict[exp_dict["name"]] = nd_tr
            self.ndata_dict[exp_dict["name"] + "_vl"] = nd_vl

            self.training["ndata"] += nd_tr
            self.validation["ndata"] += nd_vl
            self.experimental["ndata"] += nd_tr + nd_vl

        for pos_dict in self.pos_info:
            self.training["expdata"].append(pos_dict["expdata"])
            self.ndata_dict[pos_dict["name"]] = pos_dict["ndata"]
            self.training["posdatasets"].append(pos_dict["name"])

        self.ndata_dict["total_tr"] = self.training["ndata"]
        self.ndata_dict["total_vl"] = self.validation["ndata"]
        self.ndata_dict["total"] = self.experimental["ndata"]

    def _model_generation(self):
        """
        Fills the three dictionaries with the 'model' entry
            Note: before entering this function the dictionaries contain a list of inputs
                  and a list of outputs, but they are not connected.
                  This function connects inputs with outputs by injecting the PDF

        Compiles the validation and experimental models with fakes optimizers and learning rate
        as they are never trained, but this is needed by some backends in order to run evaluate on them
        """
        self.log_info("Generating the Model")

        input_list = self.input_list
        # Create the training, validation and true (data w/o replica) models:
        tr_model = MetaModel(input_list, self._pdf_injection(self.training["output"]))
        vl_model = MetaModel(input_list, self._pdf_injection(self.validation["output"]))
        true_model = MetaModel(input_list, self._pdf_injection(self.experimental["output"]))

        if self.model_file:
            # If a model file is given, load the weights from there
            # note: even though the load corresponds to the training model only, the layer_pdf is shared
            # and so it should affect all models
            tr_model.load_weights(self.model_file)

        if self.print_summary:
            tr_model.summary()

        fake_opt = "RMSprop"
        fake_lr = 0.01
        # Compile true and validation with fake lr/optimizer
        vl_model.compile(optimizer_name=fake_opt, learning_rate=fake_lr, loss=self.validation["losses"])
        true_model.compile(optimizer_name=fake_opt, learning_rate=fake_lr, loss=self.experimental["losses"])

        if self.test_dict:
            # If a test_dict was given, create and compile also the test model
            test_model = MetaModel(input_list, self._pdf_injection(self.test_dict["output"]))
            test_model.compile(optimizer_name=fake_opt, learning_rate=fake_lr, loss=self.test_dict["losses"])
            self.test_dict["model"] = test_model

        self.training["model"] = tr_model
        self.validation["model"] = vl_model
        self.experimental["model"] = true_model

    def _reset_observables(self):
        """
        Resets the 'output' and 'losses' entries of all 3 dictionaries,
        as well as the input_list
        this is necessary as these can either depend on the parametrization of the NN
        or be obliterated when the backend state is reset
        """
        self.input_list = []
        for key in ["output", "losses"]:
            self.training[key] = []
            self.validation[key] = []
            self.experimental[key] = []
            if self.test_dict:
                self.test_dict[key] = []

    def _pdf_injection(self, olist):
        """
        Takes as input a list of output layers and returns a corresponding list
        where all output layers call the pdf layer at self.pdf_layer
        """
        return [o(self.layer_pdf) for o in olist]

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
    def _generate_observables(self, params):
        """
        This functions fills the 3 dictionaries (training, validation, experimental)
        with the output layers and the loss functions
        It also fill the list of input tensors (input_list)

        Parameters accepted:
            - pos_multiplier: the multiplier to be applied to the positivity
            - pos_initial: the initial value for the positivity
        """

        # First reset the dictionaries
        self._reset_observables()
        self.log_info("Generating layers")

        # two parameters are used here which might go to the hyperopt
        pos_multiplier = params["pos_multiplier"]  # how much to increase the positivity penalty each 100 epochs
        pos_initial = params["pos_initial"]

        # Now we need to loop over all dictionaries (First exp_info, then pos_info)
        for exp_dict in self.exp_info:
            if not self.mode_hyperopt:
                self.log_info("Generating layers for experiment {0}".format(exp_dict["name"]))

            exp_layer = model_gen.observable_generator(exp_dict)

            # Save the input(s) corresponding to this experiment
            self.input_list += exp_layer["inputs"]

            if exp_dict["name"] == TESTNAME and self.test_dict:
                # If this is the test set, fill the dictionary and stop here
                self.test_dict["output"].append(exp_layer["output"])
                self.test_dict["losses"].append(exp_layer["loss"])
                continue

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
                self.log_info("Generating positivity penalty for {0}".format(pos_dict["name"]))
            pos_layer = model_gen.observable_generator(
                pos_dict, positivity_initial=pos_initial, positivity_multiplier=pos_multiplier
            )
            # The input list is still common
            self.input_list += pos_layer["inputs"]

            # The positivity all falls to the training
            self.training["output"].append(pos_layer["output_tr"])
            self.training["losses"].append(pos_layer["loss_tr"])
        # Save the positivity multiplier into the training dictionary as it will be used during training
        self.training["pos_multiplier"] = pos_multiplier

        # If there is no validation, overwrite the values of validation with the values of training

    #         if self.no_validation:
    #             self.validation['losses'] = self.training['losses']
    #             self.validation['output'] = self.training['output']

    def _generate_pdf(self, params):
        """
        Defines the internal variable layer_pdf
        this layer takes any input (x) and returns the pdf value for that x

        if the sumrule is being imposed, it also updates input_list with the
        integrator_input tensor used to calculate the sumrule

        Returns: (layers, integrator_input)
            layers: a list of layers
            integrator_input: input used to compute the  sumrule
        both are being used at the moment for reporting purposes at the end of the fit

        Parameters accepted:
            - nodes_per_layer
            - activation_per_layer
            - initializer
            - layer_type
            - dropout
        """
        self.log_info("Generating PDF layer")

        # Set the parameters of the NN
        nodes_per_layer = params["nodes_per_layer"]
        activation_per_layer = params["activation_per_layer"]
        initializer = params["initializer"]
        layer_type = params["layer_type"]
        dropout = params["dropout"]

        # Generate the NN layers
        layer_pdf, layers = model_gen.pdfNN_layer_generator(
            nodes=nodes_per_layer,
            activations=activation_per_layer,
            layer_type=layer_type,
            flav_info=self.flavinfo,
            seed=self.NNseed,
            initializer_name=initializer,
            dropout=dropout,
        )

        integrator_input = None
        if self.impose_sumrule:
            # Impose the sumrule
            # Inyect here momentum sum rule, effecively modifying layer_pdf
            layer_pdf, integrator_input = constrains.msr_impose(layers["fitbasis"], layer_pdf)
            self.input_list.append(integrator_input)

        self.layer_pdf = layer_pdf

        return layers, integrator_input

    def _model_compilation(self, params):
        """
        Compiles the model with the data given in params

        Parameters accepted:
            - learning_rate
            - optimizer
        """
        training_model = self.training["model"]
        loss_list = self.training["losses"]

        # Before compilation, set all parameters as given in the params dictionary
        learning_rate = params["learning_rate"]
        # Choose the optimizer, compile the net and print summary
        opt_name = params["optimizer"]

        training_model.compile(optimizer_name=opt_name, learning_rate=learning_rate, loss=loss_list)

    def _train_and_fit(self, validation_object, epochs):
        """
        Trains the NN for the number of epochs given using
        validation_object as the stopping criteria

        Every 100 epochs the positivitiy will be updated with
        self.training['pos_multiplier']
        """
        training_model = self.training["model"]
        exp_list = self.training["expdata"]
        pos_multiplier = self.training["pos_multiplier"]
        # Train the model for the number of epochs given
        for epoch in range(epochs):
            out = training_model.fit(x=None, y=exp_list, epochs=1, verbose=False, shuffle=False)
            passes = validation_object.monitor_chi2(out, epoch)

            if validation_object.stop_here():
                break

            if (epoch + 1) % 100 == 0:
                for pos_ds in self.training["posdatasets"]:
                    name = pos_ds
                    curr_w = training_model.get_layer(name).get_weights()
                    new_w = [curr_w[0] * pos_multiplier]
                    training_model.get_layer(name).set_weights(new_w)

        # Report a "good" training only if there was no NaNs and there was at least a point which passed positivity
        if passes and validation_object.positivity:
            return self.pass_status
        else:
            return self.failed_status

    def _hyperopt_override(self, params):
        # I love the smell of napalm in the morning
        for hyperkey in self.hyperkeys:
            item = params[hyperkey]
            if isinstance(item, dict):
                for key, value in item.items():
                    params[key] = value

    def hyperparametizable(self, params):
        """
        Wrapper around all the functions defining the fit

        This is a necessary step for hyperopt to work

        Parameters used here:
            - epochs: maximum number of iterations for the fit to run
            - stopping_patience: patience of the stopper after finding a new minimum
        """

        # Reset the internal state of the backend
        print("")
        clear_backend_state(self.debug)

        # When doing hyperopt some entries in the params dictionary can bring with them overriding arguments
        if self.mode_hyperopt:
            self.log_info("Performing hyperparameter scan")
            for key in self.hyperkeys:
                self.log_info(" > > Testing {0} = {1}".format(key, params[key]))
            self._hyperopt_override(params)

        # Fill the 3 dictionaries (training, validation, experimental) with the layers and losses
        self._generate_observables(params)

        # Generate the pdf layer
        layers, integrator_input = self._generate_pdf(params)

        # Model generation
        self._model_generation()

        # Generate the validation_object
        # this object golds statistical information about the fit and it can be used to perform stopping
        epochs = int(params["epochs"])
        stopping_patience = params["stopping_patience"]
        stopping_epochs = epochs * stopping_patience

        # If the tr/vl splitting is == 1, there will be no points in the validation so we use the training as val
        if self.no_validation:
            validation_object = Statistics.Stat_Info(
                self.training["model"],
                self.training["expdata"],
                self.ndata_dict,
                total_epochs=epochs,
                stopping_epochs=stopping_epochs,
            )
        else:
            validation_object = Statistics.Stat_Info(
                self.validation["model"],
                self.validation["expdata"],
                self.ndata_dict,
                total_epochs=epochs,
                stopping_epochs=stopping_epochs,
            )

        # Compile the training['model'] with the given parameters
        self._model_compilation(params)

        ### Training loop
        passed = self._train_and_fit(validation_object, epochs)

        # After training has completed, reload the best_model found during training
        validation_object.reload_model()

        # Compute validation loss
        validation_loss = validation_object.loss()

        # Compute experimental loss
        exp_final = self.experimental["model"].evaluate(x=None, y=self.experimental["expdata"], batch_size=1)
        try:
            experimental_loss = exp_final[0] / self.experimental["ndata"]
        except IndexError:
            experimental_loss = exp_final / self.experimental["ndata"]

        # Compute the testing loss if it was given
        if self.test_dict:
            # Generate the 'true' chi2 with the experimental model but only for models that were stopped
            target_model = self.test_dict["model"]
            target_data = self.test_dict["expdata"]
            out_final = target_model.evaluate(x=None, y=target_data, verbose=False, batch_size=1)
            try:
                testing_loss = out_final[0] / self.test_dict["ndata"]
            except IndexError:
                testing_loss = out_final / self.test_dict["ndata"]
        else:
            testing_loss = 0.0

        if self.mode_hyperopt:
            final_loss = validation_multiplier * validation_loss
            final_loss += test_multiplier * testing_loss
            final_loss /= test_multiplier + validation_multiplier
            arc_lengths = constrains.compute_arclength(layers["fitbasis"], verbose=False)
        else:
            final_loss = experimental_loss
            arc_lengths = None

        dict_out = {
            "loss": final_loss,
            "status": passed,
            "arc_lengths": arc_lengths,
            "training_loss": validation_object.tr_loss(),
            "validation_loss": validation_loss,
            "experimental_loss": experimental_loss,
            "testing_loss": testing_loss,
        }

        if self.mode_hyperopt:
            # If we are using hyperopt we don't need to output any other information
            return dict_out

        # Add to the output dictionary things that are needed by fit.py
        # to generate the output pdf, check the arc-length, gather stats, etc
        # some of them are already attributes of the class so they are redundant here
        # but I think it's good to present them explicitly
        dict_out["layer_pdf"] = self.layer_pdf
        dict_out["layers"] = layers
        dict_out["integrator_input"] = integrator_input
        dict_out["validation_object"] = validation_object
        dict_out["training"] = self.training

        return dict_out
