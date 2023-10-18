"""
    Fit action controller
"""

# Backend-independent imports
import copy
import logging

import numpy as np

import n3fit.checks
from n3fit.vpinterface import N3PDF

log = logging.getLogger(__name__)


# Action to be called by validphys
# All information defining the NN should come here in the "parameters" dict
@n3fit.checks.can_run_multiple_replicas
@n3fit.checks.check_fiatlux_pdfs_id
def performfit(
    *,
    n3fit_checks_action,  # wrapper for all checks
    replicas,  # checks specific to performfit
    replicas_nnseed_fitting_data_dict,
    posdatasets_fitting_pos_dict,
    integdatasets_fitting_integ_dict,
    theoryid,
    fiatlux,
    basis,
    fitbasis,
    sum_rules=True,
    parameters,
    replica_path,
    output_path,
    save=None,
    load=None,
    hyperscanner=None,
    hyperopt=None,
    kfold_parameters,
    tensorboard=None,
    debug=False,
    maxcores=None,
    parallel_models=False,
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
                taking as a base the fitting dictionary and the runcard's hyperscanner dictionary
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
        fiatlux: dict
            dictionary containing the params needed from LuxQED
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
        save: None, str
            model file where weights will be saved, used in conjunction with
            ``load``.
        load: None, str
            model file from which to load weights from.
        hyperscanner: dict
            dictionary containing the details of the hyperscanner
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
        parallel_models: bool
            whether to run models in parallel
    """
    from n3fit.backends import set_initial_state

    # If debug is active, the initial state will be fixed so that the run is reproducible
    set_initial_state(debug=debug, max_cores=maxcores)

    from n3fit.stopwatch import StopWatch

    stopwatch = StopWatch()

    # All potentially backend dependent imports should come inside the fit function
    # so they can eventually be set from the runcard
    from n3fit.io.writer import WriterWrapper
    from n3fit.model_trainer import ModelTrainer

    # Note: there are three possible scenarios for the loop of replicas:
    #   1.- Only one replica is being run, in this case the loop is only evaluated once
    #   2.- Many replicas being run, in this case each will have a replica_number, seed, etc
    #       and they will be fitted sequentially
    #   3.- Many replicas being run in parallel. In this case the loop will be evaluated just once
    #       but a model per replica will be generated
    #
    # In the main scenario (1) replicas_nnseed_fitting_data_dict is a list of just one element
    # case (3) is similar but the one element of replicas_nnseed_fitting_data_dict will be modified
    # to be (
    #       [list of all replica idx],
    #       one experiment with data=(replicas, ndata),
    #       [list of all NN seeds]
    #       )
    #
    n_models = len(replicas_nnseed_fitting_data_dict)
    if parallel_models and n_models != 1:
        replicas, replica_experiments, nnseeds = zip(*replicas_nnseed_fitting_data_dict)
        # Parse the experiments so that the output data contain information for all replicas
        # as the only different from replica to replica is the experimental training/validation data
        all_experiments = copy.deepcopy(replica_experiments[0])
        for i_exp in range(len(all_experiments)):
            training_data = []
            validation_data = []
            for i_rep in range(n_models):
                training_data.append(replica_experiments[i_rep][i_exp]['expdata'])
                validation_data.append(replica_experiments[i_rep][i_exp]['expdata_vl'])
            all_experiments[i_exp]['expdata'] = np.concatenate(training_data, axis=0)
            all_experiments[i_exp]['expdata_vl'] = np.concatenate(validation_data, axis=0)
        log.info(
            "Starting parallel fits from replica %d to %d",
            replicas[0],
            replicas[0] + n_models - 1,
        )
        replicas_info = [(replicas, all_experiments, nnseeds)]
    else:
        replicas_info = replicas_nnseed_fitting_data_dict

    for replica_idxs, exp_info, nnseeds in replicas_info:
        if not parallel_models or n_models == 1:
            # Cases 1 and 2 above are a special case of 3 where the replica idx and the seed should
            # be a list of just one element
            replica_idxs = [replica_idxs]
            nnseeds = [nnseeds]
            log.info("Starting replica fit %d", replica_idxs[0])

        # Generate a ModelTrainer object
        # this object holds all necessary information to train a PDF (up to the NN definition)
        the_model_trainer = ModelTrainer(
            exp_info,
            posdatasets_fitting_pos_dict,
            integdatasets_fitting_integ_dict,
            basis,
            fitbasis,
            nnseeds,
            debug=debug,
            kfold_parameters=kfold_parameters,
            max_cores=maxcores,
            model_file=load,
            sum_rules=sum_rules,
            parallel_models=n_models,
            theoryid=theoryid,
            lux_params=fiatlux,
            replicas=replica_idxs,
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

            # Note that hyperopt will not run in parallel or with more than one model _for now_
            replica_path_set = replica_path / f"replica_{replica_idxs[0]}"
            true_best = hyper_scan_wrapper(
                replica_path_set, the_model_trainer, hyperscanner, max_evals=hyperopt
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
            if parallel_models and n_models != 1:
                # If using tensorboard when running in parallel
                # dump the debugging data to the nnfit folder
                replica_path_set = replica_path
            else:
                replica_path_set = replica_path / f"replica_{replica_idxs[0]}"
            log_path = replica_path_set / "tboard"
            the_model_trainer.enable_tensorboard(log_path, weight_freq, profiling)

        #############################################################################
        # ### Fit                                                                   #
        # This function performs the actual fit, it reads all the parameters in the #
        # "parameters" dictionary, uses them to generate the NN and trains the net  #
        #############################################################################
        result = pdf_gen_and_train_function(parameters)
        stopwatch.register_ref("replica_fitted", "replica_set")

        stopping_object = result["stopping_object"]
        log.info("Stopped at epoch=%d", stopping_object.stop_epoch)

        final_time = stopwatch.stop()
        all_chi2s = the_model_trainer.evaluate(stopping_object)

        pdf_models = result["pdf_models"]
        q0 = theoryid.get_description().get("Q0")
        pdf_instances = [N3PDF(pdf_model, fit_basis=basis, Q=q0) for pdf_model in pdf_models]
        writer_wrapper = WriterWrapper(
            replica_idxs,
            pdf_instances,
            stopping_object,
            all_chi2s,
            q0**2,
            final_time,
        )
        writer_wrapper.write_data(replica_path, output_path.name, save)

        if tensorboard is not None:
            log.info("Tensorboard logging information is stored at %s", log_path)
