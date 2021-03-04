"""
    Fit action controller
"""

# Backend-independent imports
from collections import namedtuple
import logging
import numpy as np
import n3fit.checks
from n3fit.vpinterface import N3PDF

log = logging.getLogger(__name__)


# Action to be called by validphys
# All information defining the NN should come here in the "parameters" dict
@n3fit.checks.check_consistent_basis
@n3fit.checks.wrapper_check_NN
@n3fit.checks.wrapper_hyperopt
def performfit(
    *,
    genrep, # used for checks
    data, # used for checks
    replicas_nnseed_fitting_data_dict,
    posdatasets_fitting_pos_dict,
    integdatasets_fitting_integ_dict,
    theoryid,
    basis,
    fitbasis,
    sum_rules=True,
    parameters,
    replica_path,
    output_path,
    save_weights_each=None,
    save=None,
    load=None,
    hyperscan=None,
    hyperopt=None,
    kfold_parameters,
    tensorboard=None,
    debug=False,
    maxcores=None,
):
    """
        This action will (upon having read a validcard) process a full PDF fit
        for a set of replicas.

        The input to this function is provided by validphys
        and/or defined in the runcards or commandline arguments.

        This controller is provided with:
        1. Seeds generated using the replica number and the seeds defined in the runcard.
        2. Loaded datasets with replicas generated.
            2.1 Loaded positivity/integrability sets.

        The workflow of this controller is as follows:
        1. Generate a ModelTrainer object holding information to create the NN and perform a fit
            (at this point no NN object has been generated)
            1.1 (if hyperopt) generates the hyperopt scanning dictionary
                    taking as a base the fitting dictionary  and the runcard's hyperscan dictionary
        2. Pass the dictionary of parameters to ModelTrainer
                                        for the NN to be generated and the fit performed
            2.1 (if hyperopt) Loop over point 4 for `hyperopt` number of times
        3. Once the fit is finished, output the PDF grid and accompanying files

        Parameters
        ----------
            genrep: bool
                Whether or not to generate MC replicas. (Only used for checks)
            data: validphys.core.DataGroupSpec
                containing the datasets to be included in the fit. (Only used
                for checks)
            replicas_nnseed_fitting_data_dict: list[tuple]
                list with element for each replica (typically just one) to be
                fitted. Each element
                is a tuple containing the replica number, nnseed and
                ``fitted_data_dict`` containing all of the data, metadata
                for each group of datasets which is to be fitted.
            posdatasets_fitting_pos_dict: list[dict]
                list of dictionaries containing all data and metadata for each
                positivity dataset
            integdatasets_fitting_integ_dict: list[dict]
                list of dictionaries containing all data and metadata for each
                integrability dataset
            theoryid: validphys.core.TheoryIDSpec
                Theory which is used to generate theory predictions from model
                during fit. Object also contains some metadata on the theory
                settings.
            basis: list[dict]
                preprocessing information for each flavour to be fitted.
            fitbasis: str
                Valid basis which the fit is to be ran in. Available bases can
                be found in :py:mod:`validphys.pdfbases`.
            sum_rules: bool
                Whether to impose sum rules in fit. By default set to True
            parameters: dict
                Mapping containing parameters which define the network
                architecture/fitting methodology.
            replica_path: pathlib.Path
                path to the output of this run
            output_path: str
                name of the fit
            save_weights_each: None, int
                if set, save the state of the fit every ``save_weights_each``
                epochs
            save: None, str
                model file where weights will be save, used in conjunction with
                ``load``.
            load: None, str
                model file from which to load weights from.
            hyperscan: dict
                dictionary containing the details of the hyperscan
            hyperopt: int
                if given, number of hyperopt iterations to run
            kfold_parameters: None, dict
                dictionary with kfold settings used in hyperopt.
            tensorboard: None, dict
                mapping containing tensorboard settings if it is to be used. By
                default it is None and tensorboard is not enabled.
            debug: bool
                activate some debug options
            maxcores: int
                maximum number of (logical) cores that the backend should be aware of

    """
    from n3fit.backends import set_initial_state

    # If debug is active, the initial state will be fixed so that the run is reproducible
    set_initial_state(debug=debug, max_cores=maxcores)

    from n3fit.stopwatch import StopWatch

    stopwatch = StopWatch()

    # All potentially backend dependent imports should come inside the fit function
    # so they can eventually be set from the runcard
    from n3fit.ModelTrainer import ModelTrainer
    from n3fit.io.writer import WriterWrapper

    # Note: In the basic scenario we are only running for one replica and thus this loop is only
    # run once as replicas_nnseed_fitting_data_dict is a list of just one element
    stopwatch.register_times("data_loaded")
    for replica_number, exp_info, nnseed in replicas_nnseed_fitting_data_dict:
        replica_path_set = replica_path / f"replica_{replica_number}"
        log.info("Starting replica fit %s", replica_number)

        # Generate a ModelTrainer object
        # this object holds all necessary information to train a PDF (up to the NN definition)
        the_model_trainer = ModelTrainer(
            exp_info,
            posdatasets_fitting_pos_dict,
            integdatasets_fitting_integ_dict,
            basis,
            fitbasis,
            nnseed,
            debug=debug,
            save_weights_each=save_weights_each,
            kfold_parameters=kfold_parameters,
            max_cores=maxcores,
            model_file=load,
            sum_rules=sum_rules
        )

        # This is just to give a descriptive name to the fit function
        pdf_gen_and_train_function = the_model_trainer.hyperparametrizable

        # Read up the parameters of the NN from the runcard
        stopwatch.register_times("replica_set")

        ########################################################################
        # ### Hyperopt                                                         #
        # If hyperopt is active the parameters of NN will be substituted by the#
        # hyoperoptimizable variables.                                         #
        # Hyperopt will run for --hyperopt number of iterations before leaving #
        # this block                                                           #
        ########################################################################
        if hyperopt:
            from n3fit.hyper_optimization.hyper_scan import hyper_scan_wrapper

            true_best = hyper_scan_wrapper(
                replica_path_set, the_model_trainer, parameters, hyperscan, max_evals=hyperopt,
            )
            print("##################")
            print("Best model found: ")
            for k, i in true_best.items():
                print(f" {k} : {i} ")

            # In general after we do the hyperoptimization we do not care about the fit
            # so just let this die here
            break
        ####################################################################### end of hyperopt

        # Ensure hyperopt is off
        the_model_trainer.set_hyperopt(False)

        # Enable the tensorboard callback
        if tensorboard is not None:
            profiling = tensorboard.get("profiling", False)
            weight_freq = tensorboard.get("weight_freq", 0)
            log_path = replica_path_set / "tboard"
            the_model_trainer.enable_tensorboard(log_path, weight_freq, profiling)

        #############################################################################
        # ### Fit                                                                   #
        # This function performs the actual fit, it reads all the parameters in the #
        # "parameters" dictionary, uses them to generate the NN and trains the net  #
        #############################################################################
        result = pdf_gen_and_train_function(parameters)
        stopwatch.register_ref("replica_fitted", "replica_set")

        # After the fit is run we get a 'result' dictionary with the following items:
        stopping_object = result["stopping_object"]
        pdf_model = result["pdf_model"]
        true_chi2 = result["loss"]
        training = result["training"]
        log.info("Total exp chi2: %s", true_chi2)

        # Where has the stopping point happened (this is only for debugging purposes)
        print(
            """
        > > The stopping point has been at: {0} with a loss of {1}
                which it got at {2}. Stopping degree {3}
                Positivity state: {4}
                """.format(
                stopping_object.stop_epoch,
                stopping_object.vl_chi2,
                stopping_object.e_best_chi2,
                stopping_object.stopping_degree,
                stopping_object.positivity_status(),
            )
        )

        # Create a pdf instance
        pdf_instance = N3PDF(pdf_model, fit_basis=basis)

        # Generate the writer wrapper
        writer_wrapper = WriterWrapper(
            replica_number,
            pdf_instance,
            stopping_object,
            theoryid.get_description().get("Q0") ** 2,
            stopwatch.stop(),
        )

        # Now write the data down
        training_chi2, val_chi2, exp_chi2 = the_model_trainer.evaluate(stopping_object)
        writer_wrapper.write_data(
            replica_path_set, output_path.name, training_chi2, val_chi2, true_chi2
        )

        # Save the weights to some file for the given replica
        model_file = save
        if model_file:
            model_file_path = replica_path_set / model_file
            log.info(" > Saving the weights for future in %s", model_file_path)
            # Need to use "str" here because TF 2.2 has a bug for paths objects (fixed in 2.3 though)
            pdf_model.save_weights(str(model_file_path), save_format="h5")

        # If the history of weights is active then loop over it
        # rewind the state back to every step and write down the results
        for step in range(len(stopping_object.history.reloadable_history)):
            stopping_object.history.rewind(step)
            new_path = output_path / f"history_step_{step}/replica_{replica_number}"
            # We need to recompute the experimental chi2 for this point
            training_chi2, val_chi2, exp_chi2 = the_model_trainer.evaluate(stopping_object)
            writer_wrapper.write_data(new_path, output_path.name, training_chi2, val_chi2, exp_chi2)

        # So every time we want to capture output_path.name and addd a history_step_X
        # parallel to the nnfit folder

        if tensorboard is not None:
            log.info("Tensorboard logging information is stored at %s", log_path)
