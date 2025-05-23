"""
Fit action controller
"""

# Backend-independent imports
import logging

import n3fit.checks
from n3fit.vpinterface import N3PDF

log = logging.getLogger(__name__)


# Action to be called by validphys
# All information defining the NN should come here in the "parameters" dict
@n3fit.checks.check_multireplica_qed
@n3fit.checks.check_polarized_configs
def performfit(
    *,
    experiments_data,
    n3fit_checks_action,  # wrapper for all checks
    replicas,  # checks specific to performfit
    replicas_nnseed_fitting_data_dict,
    posdatasets_fitting_pos_dict,
    integdatasets_fitting_integ_dict,
    theoryid,
    fiatlux,
    basis,
    fitbasis,
    positivity_bound,
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
    double_precision=False,
    parallel_models=True,
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
        experiments_data: list[validphys.core.DataGroupSpec]
            similar to `data` but now passed as argument to `ModelTrainer`
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
        sum_rules: str
            Whether to impose sum rules in fit. By default set to True="ALL"
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
        double_precision: bool
            whether to use double precision
        parallel_models: bool
            whether to run models in parallel
    """
    from n3fit.backends import set_initial_state

    # If debug is active, the initial state will be fixed so that the run is reproducible
    set_initial_state(debug=debug, max_cores=maxcores, double_precision=double_precision)

    from n3fit.stopwatch import StopWatch

    stopwatch = StopWatch()

    # All potentially backend dependent imports should come inside the fit function
    # so they can eventually be set from the runcard
    from n3fit.io.writer import WriterWrapper
    from n3fit.model_trainer import ModelTrainer

    # Note that this can be run in sequence or in parallel
    # To do both cases in the same loop, we uniformize the replica information as:
    # - sequential: a list over replicas, each entry containing tuples of length 1
    # - parallel: a list of length 1, containing tuples over replicas
    #
    # Add inner tuples
    replicas_info = [
        ((replica,), (experiment,), (nnseed,))
        for replica, experiment, nnseed in replicas_nnseed_fitting_data_dict
    ]

    n_models = len(replicas_info)
    if parallel_models:
        # Move replicas from outer list to inner tuples
        replicas, experiments, nnseeds = [], [], []

        for replica, experiment, nnseed in replicas_info:
            replicas.extend(replica)
            experiments.extend(experiment)
            nnseeds.extend(nnseed)

        replicas_info = [(tuple(replicas), tuple(experiments), tuple(nnseeds))]
        log.info(
            "Starting parallel fits from replica %d to %d", replicas[0], replicas[0] + n_models - 1
        )
    else:
        log.info(
            "Starting sequential fits from replica %d to %d",
            replicas[0],
            replicas[0] + n_models - 1,
        )

    for replica_idxs, exp_info, nnseeds in replicas_info:
        log.info("Starting replica fit " + str(replica_idxs))

        # Generate a ModelTrainer object
        # this object holds all necessary information to train a PDF (up to the NN definition)
        the_model_trainer = ModelTrainer(
            experiments_data,
            exp_info,
            posdatasets_fitting_pos_dict,
            integdatasets_fitting_integ_dict,
            basis,
            fitbasis,
            nnseeds,
            positivity_bound,
            debug=debug,
            kfold_parameters=kfold_parameters,
            max_cores=maxcores,
            model_file=load,
            sum_rules=sum_rules,
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

        pdf_models = result["pdf_model"].split_replicas()
        q0 = theoryid.get_description().get("Q0")
        pdf_instances = [N3PDF(pdf_model, fit_basis=basis, Q=q0) for pdf_model in pdf_models]
        writer_wrapper = WriterWrapper(
            replica_idxs, pdf_instances, stopping_object, all_chi2s, theoryid, final_time
        )
        writer_wrapper.write_data(replica_path, output_path.name, save)

        if tensorboard is not None:
            log.info("Tensorboard logging information is stored at %s", log_path)
