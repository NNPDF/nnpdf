"""
    Fit action controller
"""

# Backend-independent imports
from collections import namedtuple
import sys
import logging
import os.path
import numpy as np
from reportengine.checks import make_argcheck, CheckError

import tensorflow.keras.backend as K
import matplotlib.pyplot as plt


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
        for i in range(replica_number):
            trvalseed = np.random.randint(0, pow(2, 31))

        np.random.seed(nnseed)
        for i in range(replica_number):
            nnseed = np.random.randint(0, pow(2, 31))

        np.random.seed(mcseed)
        for i in range(replica_number):
            mcseed = np.random.randint(0, pow(2, 31))
        trvalseeds.append(trvalseed)
        nnseeds.append(nnseed)
        mcseeds.append(mcseed)

    if genrep == 0:
        mcseeds = []

    Seeds = namedtuple("Seeds", ["trvlseeds", "nnseeds", "mcseeds"])
    return Seeds(trvalseeds, nnseeds, mcseeds)


@make_argcheck
def check_consistent_hyperscan_options(hyperopt, hyperscan, fitting):
    if hyperopt is not None and hyperscan is None:
        raise CheckError("A hyperscan dictionary needs to be defined when performing hyperopt")
    if hyperopt is not None and "kfold" not in hyperscan:
        raise CheckError("hyperscan::kfold key needs to be defined when performing hyperopt")
    if hyperopt is not None and fitting["genrep"]:
        raise CheckError("During hyperoptimization we cannot generate replicas (genrep=false)")


# Action to be called by valid phys
# All information defining the NN should come here in the "parameters" dict
@check_consistent_hyperscan_options
def performfit(
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
            fitting["basis"],
            nnseed,
            debug=debug,
            save_weights_each=fitting.get("save_weights_each"),
            kfold_parameters=kfold_parameters,
            max_cores=maxcores,
        )

        # Check whether we want to load weights from a file (maybe from a previous run)
        # check whether the file exists, otherwise set it to none
        # reading the data up will be done by the model_trainer
        if fitting.get("load"):
            model_file = fitting.get("loadfile")
            log.info(" > Loading the weights from previous training from %s", model_file)
            if not os.path.isfile(model_file):
                log.warning(" > Model file %s could not be found", model_file)
                model_file = None
            else:
                the_model_trainer.model_file = model_file

        # This is just to give a descriptive name to the fit function
        pdf_gen_and_train_function = the_model_trainer.hyperparametizable

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
                stopping_object.positivity_pass(),
            )
        )

        # Check the directory exist, if it doesn't, generate it
        activation_path = replica_path_set / "activationfunction"
        os.makedirs(activation_path, exist_ok=True)

        # Creates histogram of output values at the notes, to check for saturation
        # NOTE: In keras: output = activation(dot(input, kernel) + bias),
        # from https://keras.io/api/layers/core_layers/dense/
        # NOTE: calculation of a layer check explicitly and ur understanding of how we interpre
        # the relation between all paramters is correct.
        xgrid_input = np.expand_dims(np.logspace(-9, 0, num=2000).reshape(-1, 1), axis=0)
        xgrid_save = xgrid_input.reshape(-1, 1)[0::10]
        np.savetxt(activation_path / "activationxgrid.gz", xgrid_save)
        activationinputs = []
        dense_layers = []
        for num, layer in enumerate(pdf_model.layers):
            if layer.name.startswith("dense"):
                dense_layers.append(num)
        for layer in dense_layers:
            get_layer_input = K.function(
                [pdf_model.layers[1].input], [pdf_model.layers[layer].input]
            )
            layer_inputs = get_layer_input(xgrid_input)[0][0]
            kernel = pdf_model.layers[layer].kernel.numpy()
            bias = pdf_model.layers[layer].bias.numpy()
            input_to_activation_func = layer_inputs @ kernel + bias
            for node in range(pdf_model.layers[layer].output.shape[2]):
                input_single_node = input_to_activation_func[:, node]
                input_node_save = input_single_node[0::10]
                storeinfo = np.append(
                    [int(layer - min(dense_layers) + 1), int(node + 1)], input_node_save
                )
                activationinputs.append(storeinfo)
        activationinputs = np.array(activationinputs)
        np.savetxt(
            activation_path / "activationdata.gz",
            activationinputs,
            header="The first is the layer number, the second clolumn is the node, "
            + "followed by input values.",
        )

        # Creates a PDF model for export grid
        def pdf_function(export_xgrid):
            """
            Receives an array, returns the result of the PDF for said array
            """
            # First add an extra dimensions because the model works on batches
            xgrid = np.expand_dims(export_xgrid, 0)
            result = pdf_model.predict([xgrid])

            # Now remove the spurious dimension
            return np.squeeze(result, 0)

        # Intercept the gluon PDF and refit it after including extraplation data at the PDF level
        xgrid_exp = np.logspace(np.log10(4.05e-5), 0, num=100)
        xgrid_exp = np.insert(xgrid_exp, 0, 1e-6).reshape(-1, 1)  # x-value extrapolation data point
        xgrid_exp = np.expand_dims(xgrid_exp, axis=0)
        pdfs_exp = pdf_model.predict([xgrid_exp])
        flavors = [
            ("sigma", 1),
            ("g", 2),
            ("v", 3),
            ("v3", 4),
            ("v8", 5),
            ("t3", 9),
            ("t8", 10),
            ("t15", 11),
        ]
        for flavor in flavors:  # linearly indep. evol. basis pdfs
            flavor_index = flavor[1]
            pdfs_exp[0][0, flavor_index] = 2.0  # set y-value extrapolation data point
        pdf_model.compile(
            optimizer="Adadelta",
            loss="mean_squared_error",
            metrics=["accuracy"],
            target_output=pdfs_exp,
            learning_rate=1e-2,
        )
        xdict = {"input_1": xgrid_exp, "input_2": pdf_model.x_in["input_2"]}
        pdf_model.fit(x=xdict, y=None, epochs=10000, verbose=0)
        new_xgrid = np.logspace(-6, 0, num=1000).reshape(-1, 1)
        new_xgrid = np.expand_dims(new_xgrid, axis=0)
        new_pdfs = pdf_model.predict([new_xgrid])[0]
        for flavor in flavors:
            flavor_name = flavor[0]
            flavor_index = flavor[1]
            plt.figure()
            plt.plot(new_xgrid.flatten(), new_pdfs[:, flavor_index], color="b", label="fit")
            plt.plot(xgrid_exp.flatten(), pdfs_exp[0][:, 2], "r.", markersize=6, label="Data")
            plt.xscale("log")
            plt.xlabel("$x$")
            plt.ylabel(f"$x{flavor_name}(x)$")
            plt.title(f"{flavor_name} PDF")
            plt.legend(loc="best")
            plt.savefig(f"{replica_path_set}/{flavor_name}_pdf.png")
            plt.close()

        # Generate the writer wrapper
        writer_wrapper = WriterWrapper(
            replica_number,
            pdf_function,
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

    # Save the weights to some file
    if fitting.get("save"):
        model_file = fitting.get("savefile")
        log.info(" > Saving the weights for future in %s", model_file)
        training["model"].save_weights(model_file)
