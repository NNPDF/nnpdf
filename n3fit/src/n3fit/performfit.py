"""
    Fit action controller
"""

# Backend-independent imports
from collections import namedtuple
import logging
import os.path
import numpy as np
import n3fit.checks
from n3fit.vpinterface import N3PDF

log = logging.getLogger(__name__)


def initialize_seeds(replica: list, trvlseed: int, nnseed: int, mcseed: int, genrep: bool):
    """Action to initialize seeds for random number generation.
    We initialize three different seeds. The first is the seed
    used for training/validation splits, the second is used for
    initialization of the neural network's parameters and the
    final one is the monte carlo seeds for pseudodata replica
    generation.

    The generation of these seeds depend on the replica number
    in question. This dependence comes in by sampling the random
    number generator <replica number> times in the for loop.

    Parameters
    ----------
    replica: list
        A list of replica numbers to run over typically of size one
    trvlseed: int
        Seed initialization for training/validation split
    nnseed: int
        Seed for network initialization
    mcseed: int
        Seed for pseudodata replica generation
    genrep: bool

    Returns
    -------
    seeds: NamedTuple[List, List, List]
        A namedtuple of lists containing the trvalseeds, nnseeds, mcseeds
    """
    # First set the seed variables for
    # - Tr/Vl split
    # - Neural Network initialization
    # - Replica generation
    # These depend both on the seed set in the runcard and the replica number
    trvalseeds = []
    nnseeds = []
    mcseeds = []
    for replica_number in replica:
        np.random.seed(trvlseed)
        for _ in range(replica_number):
            trvalseed = np.random.randint(0, pow(2, 31))

        np.random.seed(nnseed)
        for _ in range(replica_number):
            nnseed = np.random.randint(0, pow(2, 31))

        np.random.seed(mcseed)
        for _ in range(replica_number):
            mcseed = np.random.randint(0, pow(2, 31))
        trvalseeds.append(trvalseed)
        nnseeds.append(nnseed)
        mcseeds.append(mcseed)

    if genrep == 0:
        mcseeds = []

    Seeds = namedtuple("Seeds", ["trvlseeds", "nnseeds", "mcseeds"])
    return Seeds(trvalseeds, nnseeds, mcseeds)


# Action to be called by valid phys
# All information defining the NN should come here in the "parameters" dict
@n3fit.checks.check_consistent_basis
@n3fit.checks.wrapper_check_NN
@n3fit.checks.wrapper_hyperopt
def performfit(
    fitting,
    experiments,
    t0set,
    replica,
    replica_path,
    output_path,
    theoryid,
    posdatasets,
    integdatasets=None,
    hyperscan=None,
    hyperopt=None,
    debug=False,
    maxcores=None,
):
    """
        This action will (upon having read a validcard) process a full PDF fit for a given replica.

        The input to this function is provided by validphys
        and/or defined in the runcards or commandline arguments.

        The workflow of this controller is as follows:
        1. Generates seeds using the replica number and the seeds defined in the runcard,
        2. Read up all datasets from the given experiments and create the necessary replicas
            2.1 Read up also the positivity sets
        3. Generate a ModelTrainer object holding information to create the NN and perform a fit
            (at this point no NN object has been generated)
            3.1 (if hyperopt) generates the hyperopt scanning dictionary
                    taking as a base the fitting dictionary  and the runcard's hyperscan dictionary
        4. Pass the dictionary of parameters to ModelTrainer
                                        for the NN to be generated and the fit performed
            4.1 (if hyperopt) Loop over point 4 for `hyperopt` number of times
        5. Once the fit is finished, output the PDF grid and accompanying files

        Parameters
        ----------
            fitting: dict
                dictionary with the hyperparameters of the fit
            experiments: dict
                vp list of experiments to be included in the fit
            t0set: str
                t0set name
            replica: list
                a list of replica numbers to run over (typically just one)
            replica_path: pathlib.Path
                path to the output of this run
            output_path: str
                name of the fit
            theorid: int
                theory id number
            posdatasets: list
                list of positivity datasets
            integdatasets: list
                list of integrability datasets
            hyperscan: dict
                dictionary containing the details of the hyperscan
            hyperopt: int
                if given, number of hyperopt iterations to run
            debug: bool
                activate some debug options
            maxcores: int
                maximum number of (logical) cores that the backend should be aware of
    """

    if debug:
        # If debug is active, fix the initial state this should make the run reproducible
        # (important to avoid non-deterministic multithread or hidden states)
        from n3fit.backends import set_initial_state

        set_initial_state()
    ###############

    from n3fit.stopwatch import StopWatch

    stopwatch = StopWatch()
    # All potentially backend dependent imports should come inside the fit function
    # so they can eventually be set from the runcard
    from n3fit.ModelTrainer import ModelTrainer
    from n3fit.io.writer import WriterWrapper
    from n3fit.backends import MetaModel, operations
    import n3fit.io.reader as reader

    # Loading t0set from LHAPDF
    if t0set is not None:
        t0pdfset = t0set.load_t0()
    else:
        t0pdfset = None

    trvlseed, nnseed, mcseed, genrep = [
        fitting.get(i) for i in ["trvlseed", "nnseed", "mcseed", "genrep"]
    ]

    seeds = initialize_seeds(replica, trvlseed, nnseed, mcseed, genrep)
    trvalseeds, nnseeds, mcseeds = seeds.trvlseeds, seeds.nnseeds, seeds.mcseeds

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
    if fitting.get("diagonal_basis"):
        log.info("working in diagonal basis")

    if hyperscan and hyperopt:
        kfold_parameters = hyperscan["kfold"]
        kpartitions = kfold_parameters["partitions"]
    else:
        kfold_parameters = None
        kpartitions = None

    # First loop over the experiments
    for exp in experiments:
        log.info("Loading experiment: {0}".format(exp))
        all_exp_dicts = reader.common_data_reader(
            exp,
            t0pdfset,
            replica_seeds=mcseeds,
            trval_seeds=trvalseeds,
            kpartitions=kpartitions,
            rotate_diagonal=fitting.get("diagonal_basis"),
        )
        for i, exp_dict in enumerate(all_exp_dicts):
            all_exp_infos[i].append(exp_dict)

    # Now read all the positivity datasets
    pos_info = []
    for pos_set in posdatasets:
        log.info("Loading positivity dataset %s", pos_set)
        pos_dict = reader.positivity_reader(pos_set)
        pos_info.append(pos_dict)

    if integdatasets is not None:
        integ_info = []
        for integ_set in integdatasets:
            log.info("Loading integrability dataset %s", integ_set)
            # Use the same reader as positivity observables
            integ_dict = reader.positivity_reader(integ_set)
            integ_info.append(integ_dict)
    else:
        integ_info = None

    # Note: In the basic scenario we are only running for one replica and thus this loop is only
    # run once and all_exp_infos is a list of just than one element
    stopwatch.register_times("data_loaded")
    for replica_number, exp_info, nnseed in zip(replica, all_exp_infos, nnseeds):
        replica_path_set = replica_path / f"replica_{replica_number}"
        log.info("Starting replica fit %s", replica_number)

        # Generate a ModelTrainer object
        # this object holds all necessary information to train a PDF (up to the NN definition)
        the_model_trainer = ModelTrainer(
            exp_info,
            pos_info,
            integ_info,
            fitting["basis"],
            fitting["fitbasis"],
            nnseed,
            debug=debug,
            save_weights_each=fitting.get("save_weights_each"),
            kfold_parameters=kfold_parameters,
            max_cores=maxcores,
            model_file=fitting.get("load")
        )

        # This is just to give a descriptive name to the fit function
        pdf_gen_and_train_function = the_model_trainer.hyperparametrizable

        # Read up the parameters of the NN from the runcard
        parameters = fitting.get("parameters")
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
        tboard = fitting.get("tensorboard")
        if tboard is not None:
            profiling = tboard.get("profiling", False)
            weight_freq = tboard.get("weight_freq", 0)
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
                stopping_object.epoch_of_the_stop,
                stopping_object.vl_loss,
                stopping_object.e_best_chi2,
                stopping_object.stopping_degree,
                stopping_object.positivity_status(),
            )
        )

        # Create a pdf instance
        pdf_instance = N3PDF(pdf_model)

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

    if tboard is not None:
        log.info("Tensorboard logging information is stored at %s", log_path)

    # Save the weights to some file
    model_file = fitting.get("save")
    if model_file:
        log.info(" > Saving the weights for future in %s", model_file)
        pdf_model.save_weights(model_file)
