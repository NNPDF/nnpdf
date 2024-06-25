"""
    The ModelTrainer class is the true driver around the n3fit code

    This class is initialized with all information about the NN, inputs and outputs.
    The construction of the NN and the fitting is performed at the same time when the
    hyperparametrizable method of the function is called.

    This allows to use hyperscanning libraries, that need to change the parameters of the network
    between iterations while at the same time keeping the amount of redundant calls to a minimum
"""

from collections import namedtuple
from itertools import zip_longest
import logging

import numpy as np

from n3fit import model_gen
from n3fit.backends import NN_LAYER_ALL_REPLICAS, MetaModel, callbacks, clear_backend_state
from n3fit.backends import operations as op
from n3fit.hyper_optimization.hyper_scan import HYPEROPT_STATUSES
import n3fit.hyper_optimization.penalties
import n3fit.hyper_optimization.rewards
from n3fit.hyper_optimization.rewards import HyperLoss
from n3fit.scaler import generate_scaler
from n3fit.stopping import Stopping
from n3fit.vpinterface import N3PDF, compute_phi
from validphys.core import DataGroupSpec
from validphys.photon.compute import Photon

log = logging.getLogger(__name__)

# Threshold defaults
# Any partition with a chi2 over the threshold will discard its hyperparameters
HYPER_THRESHOLD = 50.0
CHI2_THRESHOLD = 10.0
# Each how many epochs do we increase the positivitiy Lagrange Multiplier
PUSH_POSITIVITY_EACH = 100

# Each how many epochs do we increase the integrability Lagrange Multiplier
PUSH_INTEGRABILITY_EACH = 100

# See ModelTrainer::_xgrid_generation for the definition of each field and how they are generated
InputInfo = namedtuple("InputInfo", ["input", "split", "idx"])


def _pdf_injection(pdf_layers, observables, masks):
    """
    Takes as input a list of PDF layers each corresponding to one observable (also given as a list)
    And (where neded) a mask to select the output.
    Returns a list of obs(pdf).
    Note that the list of masks don't need to be the same size as the list of layers/observables
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
        experiments_data,
        exp_info,
        pos_info,
        integ_info,
        flavinfo,
        fitbasis,
        nnseeds,
        boundary_condition,
        debug=False,
        kfold_parameters=None,
        max_cores=None,
        model_file=None,
        sum_rules=None,
        theoryid=None,
        lux_params=None,
        replicas=None,
    ):
        """
        Parameters
        ----------
            experiments_data: list
                list of `validphys.core.DataGroupSpec` containing experiments
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
            nnseeds: list(int)
                the seed used to initialise the NN for each model to be passed to model_gen
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
            theoryid: validphys.core.TheoryIDSpec
                object contining info for generating the photon
            lux_params: dict
                dictionary containing the params needed from LuxQED
            replicas: list
                list with the replicas ids to be fitted
        """
        # Save all input information
        self.exp_info = list(exp_info)
        self.pos_info = [] if pos_info is None else pos_info
        self.integ_info = [] if integ_info is None else integ_info
        self.all_info = self.exp_info[0] + self.pos_info + self.integ_info
        self.boundary_condition = boundary_condition
        self.flavinfo = flavinfo
        self.fitbasis = fitbasis
        self._nn_seeds = nnseeds
        self.debug = debug
        self.all_datasets = []
        self._scaler = None
        self.theoryid = theoryid
        self.lux_params = lux_params
        self.replicas = replicas
        self.experiments_data = experiments_data

        # Initialise internal variables which define behaviour
        if debug:
            self.max_cores = 1
        else:
            self.max_cores = max_cores
        self.model_file = model_file
        self.print_summary = True
        self.mode_hyperopt = False
        self.impose_sumrule = sum_rules
        self._hyperkeys = None
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
            replica_statistic = kfold_parameters.get("replica_statistic", None)
            fold_statistic = kfold_parameters.get("fold_statistic", None)
            loss_type = kfold_parameters.get("loss_type", None)
            self._hyper_loss = HyperLoss(
                loss_type=loss_type,
                replica_statistic=replica_statistic,
                fold_statistic=fold_statistic,
                penalties_in_loss=kfold_parameters.get("penalties_in_loss", False),
            )

        # Initialize the dictionaries which contain all fitting information
        self.input_list = []
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
        self.experimental = {"output": [], "expdata": [], "ndata": 0, "model": None, "folds": []}
        self.tr_masks = []

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

    def set_hyperopt(self, hyperopt_on, keys=None):
        """Set hyperopt options on and off (mostly suppresses some printing)"""
        if keys is None:
            keys = []
        self._hyperkeys = keys
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
        for exp_dict in self.exp_info[0]:
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
                self.all_datasets.append(dataset.name)
        self.all_datasets = set(self.all_datasets)

        for pos_dict in self.pos_info:
            self.training["expdata"].append(pos_dict["expdata"])
            self.training["posdatasets"].append(pos_dict["name"])
            self.validation["expdata"].append(pos_dict["expdata"])
            self.validation["posdatasets"].append(pos_dict["name"])

        for integ_dict in self.integ_info:
            self.training["expdata"].append(integ_dict["expdata"])
            self.training["integdatasets"].append(integ_dict["name"])

    def _xgrid_generation(self):
        """
        Generates the full x-grid pertaining to the complete set of observables to be fitted.

        To first approximation, the full x-grid is a concatenation of all x-grid requested by
        all fk-tables.

        In the case of pineappl models all fktables ask for the same grid in x
        and so the input can be simplified to be a single grid for all (or most) datasets.
        However, this is not a _strict_ requirement for pineappl and was not a requirement before
        so the solution below must be kept general enough.

        Detailed implementation of the union of xgrids:
            let's assume an input [x1, x1, x1, x2, x2, x3]
            where each xi is a different grid, this will be broken into two lists:
            [x1, x2, x3] (unique grids) and [0,0,0,1,1,2] (index of the grid per dataset)
            The pdf will then be evaluated to concatenate([x1,x2,x3]) and then split (x1, x2, x3)
            Then each of the experiment, looking at the indexes, will receive one of the 3 PDFs
            The decision whether two grids (x1 and x1) are really the same is decided below

        The necessary information to redistribute the x-grid is held by a ``InputInfo`` tuple
        which is returned by this function.

        Returns
        ------
            Instance of ``InputInfo`` containing the input information necessary for the PDF model:
            - input:
                backend input layer with an array attached which is a concatenation of the unique
                inputs of the Model
                two inputs are the same if and only if they have the same shape, values and order
            - split:
                backend layer which splits the aforementioned concatenation back into the separate
                unique inputs, to be applied after the PDF is called
            - idx:
                indices of the observables to which the split PDF must be distributed
        """
        log.info("Generating the input grid")

        inputs_unique = []
        inputs_idx = []
        for igrid in self.input_list:
            for idx, arr in enumerate(inputs_unique):
                if igrid.size == arr.size and np.allclose(igrid, arr):
                    inputs_idx.append(idx)
                    break
            else:
                inputs_idx.append(len(inputs_unique))
                inputs_unique.append(igrid)

        # Concatenate the unique inputs
        input_arr = np.concatenate(inputs_unique, axis=1).T
        if self._scaler:
            # Apply feature scaling if given
            input_arr = self._scaler(input_arr)
        input_layer = op.numpy_to_input(input_arr)

        # The PDF model will be called with a concatenation of all inputs
        # now the output needs to be splitted so that each experiment takes its corresponding input
        sp_ar = [[i.shape[1] for i in inputs_unique]]
        sp_kw = {"axis": 2}
        sp_layer = op.as_layer(op.split, op_args=sp_ar, op_kwargs=sp_kw, name="pdf_split")

        return InputInfo(input_layer, sp_layer, inputs_idx)

    def _model_generation(self, xinput, pdf_model, partition, partition_idx):
        """
        Fills the three dictionaries (``training``, ``validation``, ``experimental``)
        with the ``model`` entry

        Compiles the validation and experimental models with fakes optimizers and learning rate
        as they are never trained, but this is needed by some backends
        in order to run evaluate on them.

        Before entering this function we have the input of the model
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
        so that every mode receives the same input and the output will be concatenated at the end
        the final output of the model is then (1, None, 14, n) (with n=number of parallel models).

        Parameters
        ----------
            xinput: InputInfo
                a tuple containing the input layer (with all values of x), and the information
                (in the form of a splitting layer and a list of indices) to distribute
                the results of the PDF (PDF(xgrid)) among the different observables
            pdf_model: n3fit.backend.MetaModel
                a model that produces PDF values
            partition: dict
                Only active during k-folding, information about the partition to be fitted
            partition_idx: int
                Index of the partition

        Returns
        -------
            models: dict
                dict of MetaModels for training, validation and experimental
        """
        log.info("Generating the Model")

        # For multireplica fits:
        #   The trainable part of the n3fit framework is a concatenation of all PDF models
        # We apply the Model as Layers and save for later the model (full_pdf)
        full_model_input_dict, full_pdf = pdf_model.apply_as_layer({"pdf_input": xinput.input})

        split_pdf_unique = xinput.split(full_pdf)

        # Now reorganize the uniques PDF so that each experiment receives its corresponding PDF
        split_pdf = [split_pdf_unique[i] for i in xinput.idx]
        # If we are in a kfolding partition, select which datasets are out
        training_mask = validation_mask = experimental_mask = [None]
        if partition and partition["datasets"]:
            # If we want to overfit the fold, leave the training and validation masks as [None]
            # otherwise, use the mask generated for the fold.
            # The experimental model instead is always limited to the fold
            if not partition.get("overfit", False):
                training_mask = [i[partition_idx] for i in self.training["folds"]]
                validation_mask = [i[partition_idx] for i in self.validation["folds"]]
            experimental_mask = [i[partition_idx] for i in self.experimental["folds"]]

        # Training and validation leave out the kofld dataset
        # experiment leaves out the negation
        output_tr = _pdf_injection(split_pdf, self.training["output"], training_mask)
        training = MetaModel(full_model_input_dict, output_tr)

        # Validation skips integrability and the "true" chi2 skips also positivity,
        # so we must only use the corresponding subset of PDF functions
        val_pdfs = []
        exp_pdfs = []
        for partial_pdf, obs in zip(split_pdf, self.training["output"]):
            if not obs.positivity and not obs.integrability:
                val_pdfs.append(partial_pdf)
                exp_pdfs.append(partial_pdf)
            elif not obs.integrability and obs.positivity:
                val_pdfs.append(partial_pdf)

        # We don't want to included the integrablity in the validation
        output_vl = _pdf_injection(val_pdfs, self.validation["output"], validation_mask)
        validation = MetaModel(full_model_input_dict, output_vl)

        # Or the positivity in the total chi2
        output_ex = _pdf_injection(exp_pdfs, self.experimental["output"], experimental_mask)
        experimental = MetaModel(full_model_input_dict, output_ex)

        if self.print_summary:
            training.summary()
            pdf_model = training.get_layer("PDFs")
            pdf_model.summary()
            nn_model = pdf_model.get_layer(NN_LAYER_ALL_REPLICAS)
            nn_model.summary()
            # We may have fits without sumrules imposed
            try:
                msr_model = pdf_model.get_layer("impose_msr")
                msr_model.summary()
            except ValueError:
                pass

        models = {"training": training, "validation": validation, "experimental": experimental}

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
        for key in ["output", "posmultipliers", "integmultipliers"]:
            self.training[key] = []
            self.validation[key] = []
            self.experimental[key] = []

    ############################################################################
    # # Parameterizable functions                                                #
    #                                                                          #
    # The functions defined in this block accept a 'params' dictionary which   #
    # defines the fit and the behaviours of the Neural Networks                #
    #                                                                          #
    # These are all called by the function hyperparamizable below              #
    # i.e., the most important function is hyperparametrizable, which is a     #
    # wrapper around all of these                                              #
    ############################################################################
    def _generate_observables(
        self,
        all_pos_multiplier,
        all_pos_initial,
        all_integ_multiplier,
        all_integ_initial,
        epochs,
        interpolation_points,
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

        # We need to transpose Experimental data, stacking over replicas
        experiment_data = {
            "trmask": [],
            "expdata": [],
            "expdata_vl": [],
            "invcovmat": [],
            "invcovmat_vl": [],
        }

        # Loop over datasets
        for i in range(len(self.exp_info[0])):
            # Loop over data fields
            for key, value in experiment_data.items():
                replica_data = []
                # Loop over replicas
                for replica in self.exp_info:
                    if key in ["expdata", "expdata_vl"]:
                        # Save the data with shape (ndata) instead of (1, ndata)
                        replica_data.append(replica[i][key][0])
                    else:
                        replica_data.append(replica[i][key])
                # Stack
                value.append(np.stack(replica_data))

        # Now we need to loop over all dictionaries (First exp_info, then pos_info and integ_info)
        for i, exp_dict in enumerate(self.exp_info[0]):
            if not self.mode_hyperopt:
                log.info("Generating layers for experiment %s", exp_dict["name"])

            # Stacked tr-vl mask array for all replicas for this dataset
            exp_layer = model_gen.observable_generator(
                exp_dict,
                self.boundary_condition,
                mask_array=experiment_data["trmask"][i],
                training_data=experiment_data["expdata"][i],
                validation_data=experiment_data["expdata_vl"][i],
                invcovmat_tr=experiment_data["invcovmat"][i],
                invcovmat_vl=experiment_data["invcovmat_vl"][i],
                n_replicas=len(self.replicas),
            )

            # Save the input(s) corresponding to this experiment
            self.input_list.append(exp_layer["inputs"])

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
            num_experiments = len(self.exp_info)
            replica_masks = np.stack([pos_dict["trmask"]] * num_experiments)
            training_data = np.stack([pos_dict["expdata"].flatten()] * num_experiments)

            pos_layer = model_gen.observable_generator(
                pos_dict,
                self.boundary_condition,
                positivity_initial=pos_initial,
                mask_array=replica_masks,
                training_data=training_data,
                validation_data=training_data,
                n_replicas=len(self.replicas),
            )
            # The input list is still common
            self.input_list.append(pos_layer["inputs"])

            # The positivity should be on both training and validation models
            self.training["output"].append(pos_layer["output_tr"])
            self.validation["output"].append(pos_layer["output_tr"])

            self.training["posmultipliers"].append(pos_multiplier)
            self.training["posinitials"].append(pos_initial)

        # Finally generate the integrability penalty
        for integ_dict in self.integ_info:
            if not self.mode_hyperopt:
                log.info("Generating integrability penalty for %s", integ_dict["name"])

            integrability_steps = int(epochs / PUSH_INTEGRABILITY_EACH)
            max_lambda = integ_dict["lambda"]

            integ_initial, integ_multiplier = _LM_initial_and_multiplier(
                all_integ_initial, all_integ_multiplier, max_lambda, integrability_steps
            )

            integ_layer = model_gen.observable_generator(
                integ_dict,
                self.boundary_condition,
                positivity_initial=integ_initial,
                integrability=True,
                n_replicas=len(self.replicas),
            )
            # The input list is still common
            self.input_list.append(integ_layer["inputs"])

            # The integrability all falls to the training
            self.training["output"].append(integ_layer["output_tr"])
            self.training["integmultipliers"].append(integ_multiplier)
            self.training["integinitials"].append(integ_initial)

        # Store a reference to the interpolator as self._scaler
        if interpolation_points:
            self._scaler = generate_scaler(self.input_list, interpolation_points)

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
        photons,
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
            photons: :py:class:`validphys.photon.compute.Photon`
                function to compute the photon PDF
        see model_gen.pdfNN_layer_generator for more information

        Returns
        -------
            pdf_model: MetaModel
                pdf model
        """
        log.info("Generating PDF models")
        pdf_model = model_gen.generate_pdf_model(
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
            scaler=self._scaler,
            num_replicas=len(self.replicas),
            photons=photons,
        )
        return pdf_model

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

    def _train_and_fit(self, training_model, stopping_object, epochs=100) -> bool:
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

    def _hyperopt_override(self, params):
        """Unrolls complicated hyperopt structures into very simple dictionaries"""
        # If the input contains all parameters, then that's your dictionary of hyperparameters
        hyperparameters = params.get("parameters")
        if hyperparameters is not None:
            return hyperparameters
        # Else, loop over all different keys and unroll the dictionaries within hyperparameters
        for hyperkey in self._hyperkeys:
            item = params[hyperkey]
            if isinstance(item, dict):
                params.update(item)
        return params

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

    def _filter_datagroupspec(self, datasets_partition):
        """Takes a list of all input exp datasets as :class:`validphys.core.DataGroupSpec`
        and select `DataSetSpec`s whose names are in datasets_partition.

        Parameters
        ----------
            datasets_partition: List[str]
                List with names of the datasets you want to select.

        Returns
        -------
            filtered_datagroupspec: List[validphys.core.DataGroupSpec]
                List of filtered exp datasets whose names are in datasets_partition.
        """
        filtered_datagroupspec = []

        # self.experiments_data is composed of a list of `DataGroupSpec` objects
        # These represent a group of related exp data sets
        # Loop over this list
        for datagroup in self.experiments_data:
            filtered_datasetspec = []

            # Each `DataGroupSpec` is composed by several `DataSetSpec` objects
            # `DataSetSpec` represents each exp dataset
            # Now, loop over them
            for dataset in datagroup.datasets:
                # Include `DataSetSpec`s whose names are in datasets_partition
                if dataset.name in datasets_partition:
                    filtered_datasetspec.append(dataset)

            # List of filtered experiments as `DataGroupSpec`
            filtered_datagroupspec.append(
                DataGroupSpec(name=f"{datagroup.name}_exp", datasets=filtered_datasetspec)
            )

        return filtered_datagroupspec

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

        # When doing hyperopt some entries in the params dictionary
        # can bring with them overriding arguments
        if self.mode_hyperopt:
            log.info("Performing hyperparameter scan")
            for key in self._hyperkeys:
                log.info(" > > Testing %s = %s", key, params[key])
            params = self._hyperopt_override(params)

        # Preprocess some hyperparameters
        epochs = int(params["epochs"])
        stopping_patience = params["stopping_patience"]
        stopping_epochs = int(epochs * stopping_patience)

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
            params.get("interpolation_points"),
        )
        threshold_pos = positivity_dict.get("threshold", 1e-6)
        threshold_chi2 = params.get("threshold_chi2", CHI2_THRESHOLD)

        # Initialize the chi2 dictionaries
        l_valid = []
        l_exper = []
        l_hyper = []
        # And lists to save hyperopt utilities
        pdfs_per_fold = []
        exp_models = []
        # phi evaluated over training/validation exp data
        trvl_phi_per_fold = []

        # Generate the grid in x, note this is the same for all partitions
        xinput = self._xgrid_generation()

        # Initialize all photon classes for the different replicas:
        if self.lux_params:
            photons = Photon(
                theoryid=self.theoryid, lux_params=self.lux_params, replicas=self.replicas
            )
        else:
            photons = None
        ### Training loop
        for k, partition in enumerate(self.kpartitions):
            # Each partition of the kfolding needs to have its own separate model
            # and the seed needs to be updated accordingly
            seeds = self._nn_seeds
            if k > 0:
                # generate random integers for each k-fold from the input `nnseeds`
                # we generate new seeds to avoid the integer overflow that may
                # occur when doing k*nnseeds
                rngs = [np.random.default_rng(seed=seed) for seed in seeds]
                seeds = [generator.integers(1, pow(2, 30)) * k for generator in rngs]

            # Generate the pdf model
            pdf_model = self._generate_pdf(
                params["nodes_per_layer"],
                params["activation_per_layer"],
                params["initializer"],
                params["layer_type"],
                params["dropout"],
                params.get("regularizer", None),  # regularizer optional
                params.get("regularizer_args", None),
                seeds,
                photons,
            )

            if photons:
                pdf_model.get_layer("add_photon").register_photon(xinput.input.tensor_content)

            # Model generation joins all the different observable layers
            # together with pdf model generated above
            models = self._model_generation(xinput, pdf_model, partition, k)

            # Only after model generation, apply possible weight file
            # Starting every replica with the same weights
            if self.model_file:
                log.info("Applying model file %s", self.model_file)
                pdf_model.load_identical_replicas(self.model_file)

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
                pdf_model,
                total_epochs=epochs,
                stopping_patience=stopping_epochs,
                threshold_positivity=threshold_pos,
                threshold_chi2=threshold_chi2,
            )

            # Compile each of the models with the right parameters
            for model in models.values():
                model.compile(**params["optimizer"])

            self._train_and_fit(models["training"], stopping_object, epochs=epochs)

            if self.mode_hyperopt:
                validation_loss = stopping_object.vl_chi2

                # number of active points in this fold
                # it would be nice to have a ndata_per_fold variable coming in the vp object...
                ndata = np.sum([np.count_nonzero(i[k]) for i in self.experimental["folds"]])
                # If ndata == 0 then it's the opposite, all data is in!
                if ndata == 0:
                    ndata = self.experimental["ndata"]

                # Compute experimental loss, over excluded datasets
                exp_loss_raw = models["experimental"].compute_losses()["loss"]
                experimental_loss = exp_loss_raw / ndata

                # Compute penalties per replica
                penalties = {
                    penalty.__name__: penalty(pdf_model=pdf_model, stopping_object=stopping_object)
                    for penalty in self.hyper_penalties
                }

                # Extracting the necessary data to compute phi
                # First, create a list of `validphys.core.DataGroupSpec`
                # containing only exp datasets within the held out fold
                experimental_data = self._filter_datagroupspec(partition["datasets"])

                # Compute per replica hyper losses
                hyper_loss = self._hyper_loss.compute_loss(
                    penalties=penalties,
                    experimental_loss=experimental_loss,
                    pdf_model=pdf_model,
                    experimental_data=experimental_data,
                    fold_idx=k,
                )

                # Create another list of `validphys.core.DataGroupSpec`
                # containing now exp datasets that are included in the training/validation dataset
                trvl_partitions = list(self.kpartitions)
                trvl_partitions.pop(k)
                trvl_exp_names = [
                    exp_name for item in trvl_partitions for exp_name in item['datasets']
                ]
                trvl_data = self._filter_datagroupspec(trvl_exp_names)
                # evaluate phi on training/validation exp set
                trvl_phi = compute_phi(N3PDF(pdf_model.split_replicas()), trvl_data)

                # Now save all information from this fold
                l_hyper.append(hyper_loss)
                l_valid.append(validation_loss)
                l_exper.append(experimental_loss)
                trvl_phi_per_fold.append(trvl_phi)
                pdfs_per_fold.append(pdf_model)
                exp_models.append(models["experimental"])

                if hyper_loss > self.hyper_threshold:
                    log.info(
                        "Loss above threshold (%.1f > %.1f), breaking",
                        hyper_loss,
                        self.hyper_threshold,
                    )
                    # Apply a penalty proportional to the number of folds not computed
                    pen_mul = len(self.kpartitions) - k
                    l_hyper = [i * pen_mul for i in l_hyper]
                    passed = False
                    break
                else:
                    passed = True
                    log.info("Fold %d finished, loss=%.1f, pass=%s", k + 1, hyper_loss, passed)

            # endfor

        if self.mode_hyperopt:
            # turn losses into arrays
            l_hyper = np.array(l_hyper)
            l_valid = np.array(l_valid)
            l_exper = np.array(l_exper)

            # Compute the loss over all folds for hyperopt
            final_hyper_loss = self._hyper_loss.reduce_over_folds(l_hyper)

            # Hyperopt needs a dictionary with information about the losses
            # it is possible to store arbitrary information in the trial file
            # by adding it to this dictionary
            dict_out = {
                "status": HYPEROPT_STATUSES[passed],
                "loss": final_hyper_loss,
                "validation_loss": np.average(l_valid),
                "experimental_loss": np.average(l_exper),
                "kfold_meta": {
                    "validation_losses": l_valid,
                    "trvl_losses_phi": np.array(trvl_phi_per_fold),
                    "experimental_losses": l_exper,
                    "hyper_losses": np.array(self._hyper_loss.chi2_matrix),
                    "hyper_losses_phi": np.array(self._hyper_loss.phi_vector),
                    "penalties": {
                        name: np.array(values)
                        for name, values in self._hyper_loss.penalties.items()
                    },
                },
            }
            return dict_out

        # Keep a reference to the models after training for future reporting
        self.training["model"] = models["training"]
        self.experimental["model"] = models["experimental"]
        self.validation["model"] = models["validation"]

        # In a normal run, the only information we need to output is the stopping object
        # (which contains metadata about the stopping)
        # and the pdf model (which are used to generate the PDF grids and compute arclengths)
        if not self.mode_hyperopt:
            passed = any(bool(i) for i in stopping_object.e_best_chi2)
        dict_out = {"status": passed, "stopping_object": stopping_object, "pdf_model": pdf_model}
        return dict_out
