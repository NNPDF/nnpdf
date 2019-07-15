"""
    Fit action controller
"""

# Backend-independent imports
import sys
import logging
import os.path
import numpy as np

log = logging.getLogger(__name__)

# Action to be called by validphys
# All information defining the NN should come here in the "parameters" dict
def fit(
    fitting,
    experiments,
    t0set,
    replica,
    replica_path,
    output_path,
    theoryid,
    posdatasets,
    hyperscan=None,
    hyperopt=None,
    debug=False,
    create_test_card=None,
):
    """
        This action will (upon having read a validcard) process a full PDF fit for a given replica.

        The input to this function is provided by validphys and defined in the runcards or commandline arguments.

        The workflow of this controller is as follows:
        1. Using the replica number and the seeds defined in the runcard, generates a seed for everything
        2. Read up all datasets from the given experiments and create the necessary replicas
            2.1 Read up also the positivity sets
        3. Generate a ModelTrainer object which hold all information to create the NN and perform a fit
            (at this point no NN object has been generated)
            3.1 (if hyperopt) generates the hyperopt scanning dictionary taking as a base the fitting dictionary
                              and the hyperscan dictionary defined in the runcard
        4. Pass the dictionary of parameters to ModelTrainer for the NN to be generated and the fit performed
            4.1 (if hyperopt) Loop over point 4 for `hyperopt` number of times
        5. Once the fit is finished, output the PDF grid and accompanying files

        # Arguments:
            - `fitting`: dictionary with the hyperparameters of the fit
            - `experiments`: vp list of experiments to be included in the fit
            - `t0set`: t0set name
            - `replica`: a list of replica numbers to run over (tipycally just one)
            - `replica_path`: path to the output of this run
            - `output_path`: name of the fit
            - `theorid`: theory id number
            - `posdatasets` : list of positivity datasets
            - `hyperscan`: dictionary containing the details of the hyperscan
            - `hyperopt`: if given, number of hyperopt iterations to run
            - `debug`: activate some debug options
            - `create_test_card`: read the runcard and output a new runcard with a test-set defined
    """

    # Call the script to generate a runcard with a test-set and immediately exit
    if create_test_card:
        from n3fit.hyper_optimization.create_testset import create_testset

        create_testset(experiments, runcard_file=create_test_card)
        return 0
    ###############

    if debug:
        # If debug is active, fix the initial state this should make the run reproducible
        # (important to avoid non-deterministic multithread or hidden states)
        from n3fit.backends import set_initial_state

        set_initial_state()
    ###############

    # All potentially backend dependent imports should come inside the fit function
    # so they can eventually be set from the runcard
    from n3fit.ModelTrainer import ModelTrainer
    from n3fit.io.writer import WriterWrapper 
    from n3fit.backends import MetaModel
    import n3fit.io.reader as reader
    import n3fit.msr as msr_constraints

    if hyperopt:
        import hyperopt as hyper
        import n3fit.hyper_optimization.filetrials as filetrials

        status_ok = hyper.STATUS_OK
    else:
        status_ok = "ok"

    # Loading t0set from LHAPDF
    if t0set is not None:
        t0pdfset = t0set.load_t0()

    # First set the seed variables for
    # - Tr/Vl split
    # - Neural Network initialization
    # - Replica generation
    # These depend both on the seed set in the runcard and the replica number
    trvalseeds = []
    nnseeds = []
    mcseeds = []
    for replica_number in replica:
        np.random.seed(fitting.get("trvlseed"))
        for i in range(replica_number):
            trvalseed = np.random.randint(0, pow(2, 31))

        np.random.seed(fitting.get("nnseed"))
        for i in range(replica_number):
            nnseed = np.random.randint(0, pow(2, 31))

        np.random.seed(fitting.get("mcseed"))
        for i in range(replica_number):
            mcseed = np.random.randint(0, pow(2, 31))
        trvalseeds.append(trvalseed)
        nnseeds.append(nnseed)
        mcseeds.append(mcseed)

    if fitting.get("genrep") == 0:
        mcseeds = []

    ##############################################################################
    # ### Read files
    # Loop over all the experiment and positivity datasets
    # and save them into two list of dictionaries: exp_info and pos_info
    # later on we will loop over these two lists and generate the actual layers
    # i.e., at this point there is no information about the NN
    # we are just creating dictionaries with all the necessary information
    # (experimental data, covariance matrix, replicas, etc, tr/val split)
    ##############################################################################
    all_exp_infos = [[] for _ in replica]
    # First loop over the experiments
    for exp in experiments:
        log.info("Loading experiment: {0}".format(exp))
        all_exp_dicts = reader.common_data_reader(exp, t0pdfset, replica_seeds=mcseeds, trval_seeds=trvalseeds)
        for i, exp_dict in enumerate(all_exp_dicts):
            all_exp_infos[i].append(exp_dict)

    # Now read all the positivity datasets
    pos_info = []
    for pos_set in posdatasets:
        log.info("Loading positivity dataset {0}".format(pos_set))
        pos_dict = reader.positivity_reader(pos_set)
        pos_info.append(pos_dict)

    # Note: In the basic scenario we are only running for one replica and thus this loop is only
    # run once and all_exp_infos is a list of just than one element
    for replica_number, exp_info, nnseed in zip(replica, all_exp_infos, nnseeds):
        replica_path_set = replica_path / "replica_{0}".format(replica_number)
        log.info("Starting replica fit {0}".format(replica_number))

        # Generate a ModelTrainer object
        # this object holds all necessary information to train a PDF (still no information about the NN)
        the_model_trainer = ModelTrainer(
            exp_info, pos_info, fitting["basis"], nnseed, pass_status=status_ok, debug=debug, save_history_each = fitting["save_history_each"]
        )

        # Check whether we want to load weights from a file (maybe from a previous run)
        # check whether the file exists, otherwise set it to none
        # reading the data up will be done by the model_trainer
        if fitting.get("load"):
            model_file = fitting.get("loadfile")
            log.info(" > Loading the weights from previous training from {0}".format(model_file))
            if not os.path.isfile(model_file):
                log.warning(" > Model file {0} could not be found".format(model_file))
                model_file = None
            else:
                the_model_trainer.model_file = model_file

        # This is just to give a descriptive name to the fit function
        pdf_gen_and_train_function = the_model_trainer.hyperparametizable

        # Read up the parameters of the NN from the runcard
        parameters = fitting.get("parameters")

        ########################################################################
        # ### Hyperopt                                                         #
        # If hyperopt is active the parameters of NN will be substituted by the#
        # hyoperoptimizable variables.                                         #
        # Hyperopt will run for --hyperopt number of iterations before leaving #
        # this block                                                           #
        ########################################################################
        if hyperopt:
            from n3fit.hyper_optimization.HyperScanner import HyperScanner

            the_scanner = HyperScanner(parameters)

            stopping_options = hyperscan.get("stopping")
            positivity_options = hyperscan.get("positivity")
            optimizer_options = hyperscan.get("optimizer")
            architecture_options = hyperscan.get("architecture")

            # Enable scanner for certain parameters
            the_scanner.stopping(**stopping_options)
            the_scanner.positivity(**positivity_options)
            the_scanner.optimizer(**optimizer_options)
            # Enable scanner for specific architectures
            the_scanner.NN_architecture(**architecture_options)

            # Tell the trainer we are doing hyperopt
            the_model_trainer.set_hyperopt(on=True, keys=the_scanner.hyper_keys)

            # Generate Trials object
            trials = filetrials.FileTrials(replica_path_set, log=log, parameters=parameters)

            # Perform the scan
            try:
                best = hyper.fmin(
                    fn=pdf_gen_and_train_function,
                    space=the_scanner.dict(),
                    algo=hyper.tpe.suggest,
                    max_evals=hyperopt,
                    trials=trials,
                )
            except ValueError as e:
                print("Error from hyperopt because no best model was found")
                print("@fit.py, setting the best trial to empty dict")
                print("Exception: {0}".format(e))
                sys.exit(0)

            # Now update the parameters with the ones found by the scan
            true_best = hyper.space_eval(parameters, best)
            the_scanner.update_dict(true_best)
            parameters = the_scanner.dict()
            print("##################")
            print("Best model found: ")
            for k, i in true_best.items():
                print(" {0} : {1} ".format(k, i))
        ####################################################################### end of hyperopt

        # Ensure hyperopt is off
        the_model_trainer.set_hyperopt(on=False)

        #############################################################################
        # ### Fit                                                                   #
        # This function performs the actual fit, it reads all the parameters in the #
        # "parameters" dictionary, uses them to generate the NN and trains the net  #
        #############################################################################
        result = pdf_gen_and_train_function(parameters)

        # After the fit is run we get a 'result' dictionary with the following items:
        validation_object = result["validation_object"]
        layer_pdf = result["layer_pdf"]
        layers = result["layers"]
        integrator_input = result["integrator_input"]
        true_chi2 = result["loss"]
        training = result["training"]
        log.info("Total exp chi2: {0}".format(true_chi2))

        # Where has the stopping point happened (this is only for debugging purposes)
        print(
            """
        > > The stopping point has been at: {0} with a loss of {1}
                which it got at {2}. Stopping degree {3}
                Positivity state: {4}
                """.format(
                validation_object.epoch_of_the_stop,
                validation_object.loss,
                validation_object.e_best_chi2,
                validation_object.stopping_degree,
                validation_object.positivity_pass(),
            )
        )

        # Creates a PDF model for export grid
        def pdf_function(export_xgrid):
            """
            Receives an array, returns the result of the PDF for said array
            """
            modelito = MetaModel([integrator_input], [], extra_tensors=[(export_xgrid, layer_pdf)])
            result = modelito.predict(x=None, steps=1)
            return result

        # Generate the writer wrapper
        writer_wrapper = WriterWrapper(pdf_function, validation_object, layers['fitbasis'], theoryid.get_description().get("Q0") ** 2)

        # Now write the data down
        writer_wrapper.write_data(replica_number, replica_path_set, output_path.stem, true_chi2)

        # If the history is active, loop over it writing down the data to different paths
        for i in validation_object.history_loop():
            # Each step of the loop reloads a different point in history
            new_path = output_path.stem + "/history_step_{0}/replica_{1}".format(i, replica_number)
            # We need to recompute the experimental chi2 for this point
            exp_chi2 = result["experimental"]["model"].evaluate()[0] / result["experimental"]["ndata"]
            writer_wrapper.write_data(replica_number, new_path, output_path.stem, exp_chi2)
            
        # So every time we want to capture output_path.stem and addd a history_step_X 
        # parallel to the nnfit folder


    # Save the weights to some file
    if fitting.get("save"):
        model_file = fitting.get("savefile")
        log.info(" > Saving the weights for future in {0}".format(model_file))
        training["model"].save_weights(model_file)

    # Plot the validation and the training losses
    if fitting.get("plot"):
        validation_object.plot()

    # print out the integration of the sum rule in case we want to check it's not broken
    # msr_constraints.check_integration(layer_pdf, integrator_input)
