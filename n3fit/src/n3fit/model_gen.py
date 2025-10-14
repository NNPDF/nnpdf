"""
Library of functions which generate the models used by n3fit to determine PDF.

It contains functions to generate:

1) Observables
    The main function is ``observable_generator`` which takes the input theory
    and generates the path from the PDF result to the computation of the
    training and validation losses / chi2

2) PDFs
    The main function is ``generate_pdf_model``, which takes a list of settings
    defining the replica-dependent architecture of each of the models that form
    the ensemble as well as ensemble-wide options such as the flavour basis,
    sum rule definition or theoretical settings, and generates a PDF model
    which takes an array of (x) as input and outputs the value of the PDF
    for each replica, for each x for each flavour.
"""

from dataclasses import asdict, dataclass, field
from typing import Callable

import numpy as np
import scipy.linalg as la

from n3fit.backends import (
    NN_LAYER_ALL_REPLICAS,
    NN_PREFIX,
    PREPROCESSING_LAYER_ALL_REPLICAS,
    Input,
    Lambda,
    MetaLayer,
    MetaModel,
    base_layer_selector,
)
from n3fit.backends import operations as op
from n3fit.layers import (
    DIS,
    DY,
    AddPhoton,
    FkRotation,
    FlavourToEvolution,
    Mask,
    ObsRotation,
    Preprocessing,
    losses,
)
from n3fit.layers.observable import is_unique
from n3fit.msr import generate_msr_model_and_grid
from validphys.photon.compute import Photon

from n3fit.backends import regularizer_selector  # isort: skip isort and black don't agree


@dataclass
class ObservableWrapper:
    """Wraps many observables into an experimental layer once the PDF model is prepared
    It can take normal datasets or Lagrange-multiplier-like datasets
    (such as positivity or integrability)
    """

    # IDEALLY:
    # In principle this is something that could be automatically provided by validphyts
    # __but__ it requires backend (i.e., tensorflow) information
    # but maybe it can be constructed in such a way that the backend is lazyly imported
    # and make this part of the experiment spec

    name: str
    observables: list
    trvl_mask_layer: Mask
    dataset_xsizes: list
    sqrtcov: np.array = None
    multiplier: float = 1.0
    integrability: bool = False
    positivity: bool = False
    data: np.array = None
    rotation: ObsRotation = None  # only used for diagonal covmat

    def _generate_loss(self, mask=None):
        """Generates the corresponding loss function depending on the values the wrapper
        was initialized with"""
        if self.sqrtcov is not None:
            sqrtcov = self.sqrtcov
            loss = losses.LossInvcovmat(self.data, mask, sqrtcov=sqrtcov, name=self.name)
        elif self.positivity:
            loss = losses.LossPositivity(name=self.name, c=self.multiplier)
        elif self.integrability:
            loss = losses.LossIntegrability(name=self.name, c=self.multiplier)
        return loss

    def _generate_experimental_layer(self, pdf):
        """Generate the experimental layer by feeding to each observable its PDF.
        In the most general case, each observable might need a PDF evaluated on a different xgrid,
        the input PDF is evaluated in all points that the experiment needs and needs to be split
        """
        if len(self.dataset_xsizes) > 1:
            splitting_layer = op.tensor_splitter(
                pdf.shape, self.dataset_xsizes, axis=2, name=f"{self.name}_split"
            )
            sp_pdf = splitting_layer(pdf)
            output_layers = [obs(p) for obs, p in zip(self.observables, sp_pdf)]
        else:
            output_layers = [obs(pdf) for obs in self.observables]

        # Finally concatenate all observables (so that experiments are one single entity)
        ret = op.concatenate(output_layers, axis=-1)

        # rotate the predictions to the diagonal basis if rotation is provided.
        if self.rotation is not None:
            ret = self.rotation(ret)

        if self.trvl_mask_layer is not None:
            ret = self.trvl_mask_layer(ret)

        return ret

    def __call__(self, pdf_layer, mask=None):
        loss_f = self._generate_loss(mask)
        experiment_prediction = self._generate_experimental_layer(pdf_layer)
        return loss_f(experiment_prediction)


def observable_generator(
    spec_dict,
    boundary_condition=None,
    training_mask_array=None,
    validation_mask_array=None,
    training_data=None,
    validation_data=None,
    covmat_tr=None,
    covmat_vl=None,
    positivity_initial=1.0,
    integrability=False,
    n_replicas=1,
):  # pylint: disable=too-many-locals
    """
    This function generates the observable models for each experiment.
    These are models which takes as input a PDF tensor (1 x size_of_xgrid x flavours) and outputs
    the result of the observable for each contained dataset (n_points,).

    In summary the model has the following structure:
        Observable layers, corresponding to commondata datasets
        and made of any number of fktables (and an operation on them).

    An observable contains an fktable, which is loaded by the convolution layer
    (be it hadronic or DIS) and a inv covmat which loaded by the loss.

    This function also outputs three "output objects" (which are functions that generate layers)
    that use the training and validation mask to create a training_output, validation_output
    and experimental_output

    If the dataset is a positivity dataset acts in consequence.

    The output is a dictionary (`layer_info`), each one of the three output functions
    have a signature:

        `def out_tr(pdf_layer, dataset_out=None)`

    The `pdf_layer` must be a layer of shape (1, size_of_xgrid, flavours)
    `datasets_out` is the list of dataset to be masked to 0 when generating the layer

    Parameters
    ----------
        spec_dict: dict
            a dictionary-like object containing the information of the experiment
        boundary_condition: dict
            dictionary containing the instance of the a PDF set to be used as a
            Boundary Condition.
        training_mask_array: np.ndarray
            training mask per replica
        validation_mask_array: np.ndarray
            validation mask per replica, when not given ~training_mask_array will be used
            while in general the validation is a negation of the training, in special cases
            such as 1-point datasets, these are accepted by both masks and then removed by the loss
        n_replicas: int
            number of replicas fitted simultaneously
        positivity_initial: float
            set the positivity lagrange multiplier for epoch 1
        integrability: bool
            switch on/off the integrability constraints

    Returns
    ------
        layer_info: dict
            a dictionary with:
            - `inputs`: input layer
            - `output`: output layer (unmasked)
            - `output_tr`: output layer (training)
            - `output_vl`: output layer (validation)
            - `experiment_xsize`: int (size of the output array)
    """
    spec_name = spec_dict["name"]
    dataset_xsizes = []
    model_inputs = []
    model_observables = []
    # The first step is to compute the observable for each of the datasets
    for dataset in spec_dict["datasets"]:
        # Get the generic information of the dataset
        dataset_name = dataset.name

        # Look at what kind of layer do we need for this dataset
        if dataset.hadronic:
            Obs_Layer = DY
        else:
            Obs_Layer = DIS

        # Set the operation (if any) to be applied to the fktables of this dataset
        operation_name = dataset.operation

        # Now generate the observable layer, which takes the following information:
        # operation name
        # dataset name
        # list of validphys.coredata.FKTableData objects
        #   these will then be used to check how many different pdf inputs are needed
        #   (and convolutions if given the case)
        obs_layer = Obs_Layer(
            dataset.fktables_data,
            dataset.fktables(),
            dataset_name,
            boundary_condition,
            operation_name,
            n_replicas=n_replicas,
            name=f"dat_{dataset_name}",
        )

        # If the observable layer found that all input grids are equal, the splitting will be None
        # otherwise the different xgrids need to be stored separately
        # Note: for pineappl grids, obs_layer_tr.splitting should always be None
        if obs_layer.splitting is None:
            xgrid = dataset.fktables_data[0].xgrid
            model_inputs.append(xgrid)
            dataset_xsizes.append(len(xgrid))
        else:
            xgrids = [i.xgrid for i in dataset.fktables_data]
            model_inputs += xgrids
            dataset_xsizes.append(sum([len(i) for i in xgrids]))

        model_observables.append(obs_layer)

    # Check whether all xgrids of all observables in this experiment are equal
    # if so, simplify the model input
    if is_unique(model_inputs):
        model_inputs = model_inputs[0:1]
        dataset_xsizes = dataset_xsizes[0:1]

    # Reshape all inputs arrays to be (1, nx)
    model_inputs = np.concatenate(model_inputs).reshape(1, -1)

    # Make the mask layers...
    if training_mask_array is None:
        tr_mask_layer = None
        if validation_mask_array is None:
            vl_mask_layer = None
        else:
            vl_mask_layer = Mask(validation_mask_array, name=f"vlmask_{spec_name}")
    else:
        tr_mask_layer = Mask(training_mask_array, name=f"trmask_{spec_name}")
        if validation_mask_array is None:
            vl_mask_layer = Mask(~training_mask_array, name=f"vlmask_{spec_name}")
        else:
            vl_mask_layer = Mask(validation_mask_array, name=f"vlmask_{spec_name}")

    # get rotation matrix to diagonal basis
    if spec_dict.get("data_transformation") is not None:
        obsrot = ObsRotation(spec_dict.get("data_transformation"))
    else:
        obsrot = None

    if spec_dict["positivity"]:
        out_positivity = ObservableWrapper(
            spec_name,
            model_observables,
            tr_mask_layer,
            dataset_xsizes,
            multiplier=positivity_initial,
            positivity=not integrability,
            integrability=integrability,
        )

        layer_info = {
            "inputs": model_inputs,
            "output_tr": out_positivity,
            "experiment_xsize": sum(dataset_xsizes),
        }
        # For positivity we end here
        return layer_info

    out_tr = ObservableWrapper(
        spec_name,
        model_observables,
        tr_mask_layer,
        dataset_xsizes,
        sqrtcov=[la.cholesky(cov, lower=True) for cov in covmat_tr],
        data=training_data,
        rotation=obsrot,
    )

    out_vl = ObservableWrapper(
        f"{spec_name}_val",
        model_observables,
        vl_mask_layer,
        dataset_xsizes,
        sqrtcov=[la.cholesky(cov, lower=True) for cov in covmat_vl],
        data=validation_data,
        rotation=obsrot,
    )

    # experimental data has already been rotated if diagonal basis is requested
    out_exp = ObservableWrapper(
        f"{spec_name}_exp",
        model_observables,
        None,
        dataset_xsizes,
        sqrtcov=la.cholesky(spec_dict["covmat"], lower=True),
        data=spec_dict["expdata_true"],
        rotation=None,
    )

    layer_info = {
        "inputs": model_inputs,
        "output": out_exp,
        "output_tr": out_tr,
        "output_vl": out_vl,
        "experiment_xsize": sum(dataset_xsizes),
    }
    return layer_info


@dataclass
class ReplicaSettings:
    """Dataclass which holds all necessary replica-dependent information of a PDF.

    Parameters
    ----------
        seed: int
            seed for the initialization of the neural network
        nodes: list[int]
            nodes of each of the layers, starting at the first hidden layer
        activations: list[str]
            list of activation functions, should be of equal length as nodes
        architecture: str
            select the architecture of the neural network used for the replica,
            e.g. ``dense`` or ``dense_per_flavour``
        initializer: str
            initializer to be used for this replica
        dropout: float
            rate of dropout for each layer
        regularizer: str
            name of the regularizer to use for this replica (if any)
        regularizer_args: dict
            options to pass down to the regularizer (if any)
    """

    seed: int
    nodes: list[int]
    activations: list[str]
    architecture: str = "dense"
    initializer: str = "glorot_normal"
    dropout_rate: float = 0.0
    regularizer: str = None
    regularizer_args: dict = field(default_factory=dict)

    def __post_init__(self):
        """Apply checks to the input, and expand hyperopt callables"""
        # Expansions
        if callable(self.activations):
            # Hyperopt might pass down a function to generate the list of activations
            # depending on the number of layers
            self.activations = self.activations(len(self.nodes))

        if self.regularizer_args is None:
            self.regularizer_args = dict()

        # Checks
        if len(self.nodes) != len(self.activations):
            raise ValueError(
                f"nodes and activations do not match ({self.nodes} vs {self.activations}"
            )
        if self.regularizer_args and self.regularizer is None:
            raise ValueError(
                "Regularizer arguments have been provided but no regularizer is selected"
            )


def generate_pdf_model(
    replicas_settings: list[ReplicaSettings],
    flav_info: dict = None,
    fitbasis: str = "NN31IC",
    out: int = 14,
    impose_sumrule: str = None,
    scaler: Callable = None,
    photons: Photon = None,
):
    """
    Generation of the full PDF model which will be used to determine the full PDF.
    The full PDF model can have any number of replicas, which can be trained in parallel,
    the limitations of the determination means that there are certain traits that all replicas
    must share, while others are fre per-PDF.

    In its most general form, the output of this function is a :py:class:`n3fit.backend.MetaModel`
    with the following architecture:

        <input layer>
            in the standard PDF fit this includes only the (x) grid of the NN

        [ list of a separate architecture per replica ]
            which can be, but is not necessary, equal for all replicas

        [ <preprocessing factors> ]
            postprocessing of the network output by a variation x^{alpha}*(1-x)^{beta}

        <normalization>
            physical sum rules, requires an integral over the PDF

        <rotation to FK-basis>
            regardless of the physical basis in which the PDF and preprocessing factors are applied
            the output is rotated to the 14-flavour general basis used in FkTables following
            PineaAPPL's convention

        [<output layer>]
            14 flavours per value of x per replica
            note that, depending on the fit basis (and fitting scale)
            the output of the PDF will contain repeated values


    This function defines how the PDFs will be generated.
    In the case of identical PDF models (``identical_models = True``, default) the same
    settings will be used for all replicas.
    Otherwise, the sampling routines will be used.


    Parameters:
    -----------
        replica_settings: list[ReplicaSettings]
            list of ReplicaSettings objects which must contain the following information
                nodes: list(int)
                    list of the number of nodes per layer of the PDF NN
                activation: list
                    list of activation functions to apply to each layer
                initializer_name: str
                    selects the initializer of the weights of the NN. Default: glorot_normal
                layer_type: str
                    selects the type of architecture of the NN. Default: dense
                seed: int
                    the initialization seed for the NN
                dropout: float
                    rate of dropout layer by layer
                regularizer: str
                    name of the regularizer to use for the NN
                regularizer_args: dict
                    options to pass down to the regularizer (if any)
                flav_info: dict
                    dictionary containing the information about each PDF (basis dictionary in the runcard)
                    to be used by Preprocessing
                fitbasis: str
                    fitbasis used during the fit. Default: NN31IC
        out: int
            number of output flavours of the model (default 14)
        impose_sumrule: str
            whether to impose sumrules on the output pdf and which one to impose (All, MSR, VSR, TSR)
        scaler: callable
            Function to apply to the input. If given the input to the model
            will be a (1, None, 2) tensor where dim [:,:,0] is scaled
            When None, instead turn the x point into a (x, log(x)) pair
        photons: :py:class:`validphys.photon.compute.Photon`
            If given, gives the AddPhoton layer a function to compute a photon which will be added at the
            index 0 of the 14-size FK basis
            This same function will also be used to compute the MSR component for the photon

    Returns
    -------
        pdf_model: MetaModel
            pdf model, with `single_replica_generator` attached as an attribute
    """
    shared_config = {
        "flav_info": flav_info,
        "fitbasis": fitbasis,
        "output_size": out,
        "impose_sumrule": impose_sumrule,
        "scaler": scaler,
        "photons": photons,
    }

    pdf_model = _pdfNN_layer_generator(replicas_settings, **shared_config)

    # Note that the photons are passed unchanged to the single replica generator
    # computing the photon requires running fiatlux which takes 30' per replica
    # and so at the moment parallel photons are disabled with a check in checks.py
    # In order to enable it `single_replica_generator` must take the index of the replica
    # to select the appropiate photon as all of them will be computed and fixed before the fit

    def single_replica_generator(replica_idx=0):
        """Generate one single replica from the entire batch.
        The select index is relative to the batch, not the entire PDF determination.

        This function is necessary to separate all the different models after training.
        """
        settings = replicas_settings[replica_idx]
        # TODO:
        # In principle we want to recover the initial replica exactly,
        # however, for the regression tests to pass
        # _in the polarized case and only in the polarized case_ this line is necessary
        # it most likely has to do with numerical precision, but panicking might be in order
        settings.seed = 0
        return _pdfNN_layer_generator([settings], **shared_config, replica_axis=False)

    pdf_model.single_replica_generator = single_replica_generator

    return pdf_model


def _pdfNN_layer_generator(
    replicas_settings: list[ReplicaSettings],
    flav_info: dict = None,
    fitbasis: str = "NN31IC",
    output_size: int = 14,
    impose_sumrule: str = None,
    scaler: Callable = None,
    photons: Photon = None,
    replica_axis: bool = True,
):
    """
    Generates the PDF model which takes as input a point in x (from 0 to 1)
    and outputs a basis of 14 PDFs.
    It generates the preprocessing of the x into a set (x, log(x)),
    the arbitrary NN to fit the form of the PDF
    and the preprocessing factor.

    The funtional form of the output of this function is of:

        f_{i}(x) = R_{ji} N(x)_{j} * x^{1-alpha_{j}} * (1-x)^{beta_{j}}

    Where i goes from 1 to 14 while j goes from 1 to the size of the basis. R_{ji}
    is the rotation from the fitting basis to the physical basis needed for the
    convolution with the fktables.

    `layer_type` defines the architecture of the Neural Network, currently
    the following two options are implemented:
        - `dense`
            it will generate a densely connected networks where all nodes of layer n
            are connected with all nodes of layer n+1 (and n-1), the output layer is
            of the size of the last item in the `nodes` list (usually 8)
        - `dense_per_flavour`
            similar to `dense` but the nodes are disconnected in a per-flavour basis
            This means at the end of the PDF (and before the preprocessing)
            we have 8 (or whatever number) nodes just as before, but from the input
            to the output the nodes are disconnected

    This is a complicated function comprised of several sections:

    1. As initialization a number of checks are carried outs to ensure the number of layers
    and the number of activation functions match up.

    2. The neural network NN(x) is constructed depending on the requested layer_type

    3. Break the input into (x, log(x)) and passes down to the NN.
    The output of the NN in the previous step is a list of independent layers.
    A function is constructed that joins all those layers. The function takes a
    tensor as the input and applies all layers for NN in order.

    4. Create a preprocessing factor layer (that takes as input the same tensor x as the NN)
    and multiply it to the NN. We have now:
        N(x)_{j} * x^{1-alpha_{j}} * (1-x)^{beta_{j}}

    5. Create a rotation layer (that takes as input the previous step) and applies it
    so that we get f(x)_{i}

    Finally we output the final answer as well as the list of all generating functions
    in the model for easy usage within `n3fit`.

    Example
    -------

    >>> import numpy as np
    >>> from n3fit.vpinterface import N3PDF
    >>> from n3fit.model_gen import _pdfNN_layer_generator, ReplicaSettings
    >>> from validphys.pdfgrids import xplotting_grid
    >>> rp = [ReplicaSettings(nodes = [8], activations=["linear"], seed=i) for i in [1,2]]
    >>> fake_fl = [{'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'g', 's', 'sbar']]
    >>> fake_x = np.linspace(1e-3,0.8,3).reshape(1,-1,1)
    >>> pdf_model = _pdfNN_layer_generator(rp, flav_info=fake_fl, fitbasis='FLAVOUR', impose_sumrule=False)
    >>> pdf_model(fake_x).shape
    TensorShape([1, 2, 3, 14])

    # 1 batch, 2 replicas, 3 x points, 14 flavours


    Parameters
    ----------
        replicas_settings: list(:py:class:`ReplicaSettings`)
            list of ``ReplicaSettings`` objects holding the settings of each of the replicas
        flav_info: dict
            dictionary containing the information about each PDF (basis dictionary in the runcard)
            to be used by Preprocessing
        fitbasis: str
            fitbasis used during the fit. Default: NN31IC
        output_size: int
            number of output flavours of the model (default 14)
        impose_sumrule: str
            whether to impose sumrules on the output pdf and which one to impose (All, MSR, VSR, TSR)
        scaler: callable
            Function to apply to the input. If given the input to the model
            will be a (1, None, 2) tensor where dim [:,:,0] is scaled
            When None, instead turn the x point into a (x, log(x)) pair
        photon: :py:class:`validphys.photon.compute.Photon`
            If given, gives the AddPhoton layer a function to compute a photon which will be added at the
            index 0 of the 14-size FK basis
            This same function will also be used to compute the MSR component for the photon
        replica_axis: bool
            Whether there is an explicit replica axis. True even with one replica for the main model,
            used internally for the single replica models.

    Returns
    -------
       pdf_model: n3fit.backends.MetaModel
            a model f(x) = y where x is a tensor (1, xgrid, 1) and y a tensor (1, replicas, xgrid, out)
    """
    all_seed = [i.seed for i in replicas_settings]
    num_replicas = len(replicas_settings)

    if impose_sumrule is None:
        impose_sumrule = "All"

    ## Process the input data (x grid)
    # There a currently two options:
    # 1. Append log(x) to the input
    # 2. Scale the input
    do_nothing = lambda x: x
    model_input = {}

    if scaler is None:  # add log(x)
        use_feature_scaling = subtract_one = False
        # The PDF itself receives only x
        pdf_input_dimensions = 1
        # But the NN will see (x, log(x))
        nn_input_dimensions = 2

        pdf_input = Input(shape=(None, pdf_input_dimensions), batch_size=1, name="pdf_input")
        process_input = Lambda(lambda x: op.concatenate([x, op.op_log(x)], axis=-1), name="x_logx")
        extract_original = do_nothing
        extract_nn_input = do_nothing
    else:
        use_feature_scaling = subtract_one = True
        # The NN will only receive x
        nn_input_dimensions = 1
        # But the PDF itself will receive both (x, scaler(x))
        pdf_input_dimensions = 2

        pdf_input = Input(shape=(None, pdf_input_dimensions), batch_size=1, name="scaledx_x")
        process_input = do_nothing
        extract_nn_input = Lambda(lambda x: op.op_gather_keep_dims(x, 0, axis=-1), name="x_scaled")
        extract_original = Lambda(lambda x: op.op_gather_keep_dims(x, 1, axis=-1), name="pdf_input")

    if subtract_one:
        # TODO: make it its own option, even though now it only activates in the scaler if above
        input_x_eq_1 = [1.0]
        if use_feature_scaling:
            input_x_eq_1 = scaler(input_x_eq_1)[0]
        # the layer that subtracts 1 from the NN output
        subtract_one_layer = Lambda(op.op_subtract, name="subtract_one")
        layer_x_eq_1 = op.numpy_to_input(np.array(input_x_eq_1).reshape(1, 1), name="x_eq_1")
        model_input["layer_x_eq_1"] = layer_x_eq_1

    model_input["pdf_input"] = pdf_input

    ## Create the actual NeuralNetwork PDF
    # loop over the settings for all replicas and generate a list of NN per replica
    # which will be then stack together and built into a single (input -> output) MetaModel
    # all PDFs _must_ share the same input layer
    x_input = Input(shape=(None, nn_input_dimensions), batch_size=1, name="NN_input")

    list_of_nn_pdfs = []
    for i, replica_settings in enumerate(replicas_settings):
        rep_pdf = _generate_nn(x_input, i, **asdict(replica_settings))
        # And build them all with the same input layer
        list_of_nn_pdfs.append(rep_pdf(x_input))

    # Stack all replicas together as one single object
    nn_pdfs = Lambda(lambda nns: op.stack(nns, axis=1), name="stack_replicas")(list_of_nn_pdfs)
    nn_replicas = MetaModel({'NN_input': x_input}, nn_pdfs, name=NN_LAYER_ALL_REPLICAS)

    ## Preprocessing factors:
    # the layer that multiplies the NN output by the preprocessing factor
    # This includes
    #       - x^{a}(1-x)^{b}
    #       - NN(x) - N(1.0)
    apply_preprocessing_factor = Lambda(op.op_multiply, name="prefactor_times_NN")

    compute_preprocessing_factor = Preprocessing(
        flav_info=flav_info,
        name=PREPROCESSING_LAYER_ALL_REPLICAS,
        replica_seeds=all_seed,
        large_x=not subtract_one,
    )

    # The NN subtracted by NN(1), if applicable, otherwise do nothing
    def nn_subtracted(x):
        NNs_x = nn_replicas(x)

        if subtract_one:
            x_eq_1_processed = process_input(layer_x_eq_1)
            NNs_x_1 = nn_replicas(x_eq_1_processed)
            NNs_x = subtract_one_layer([NNs_x, NNs_x_1])

        return NNs_x

    ## Unnormalized PDF
    #   updf_r(x) = FkRotation( NN_r(input(x)) * preprocessing_layer_r(x) )
    #       with _r: replica index
    #       input: whatever processing is applied to the input
    # The preprocessing_layer and weights is specific to each replica
    # The final PDF will be in the 14 flavours evolution basis used in the FkTables

    # Basis rotation
    basis_rotation = FlavourToEvolution(
        flav_info=flav_info, fitbasis=fitbasis, name="pdf_evolution_basis"
    )

    # Evolution layer
    layer_evln = FkRotation(output_dim=output_size, name="pdf_FK_basis")

    def compute_unnormalized_pdf(x):
        # Preprocess the input grid
        x_nn_input = extract_nn_input(x)
        x_processed = process_input(x_nn_input)
        x_original = extract_original(x)

        # Compute the neural network output
        NNs_x = nn_subtracted(x_processed)

        # Compute the preprocessing factor
        preprocessing_factors_x = compute_preprocessing_factor(x_original)

        # Apply the preprocessing factor
        pref_NNs_x = apply_preprocessing_factor([preprocessing_factors_x, NNs_x])

        # Transform to FK basis, this is the full evolution basis
        # Rotate to the 9f evolution basis first before expanding up to 14f
        # TODO: make these two steps into one
        if not basis_rotation.is_identity():
            pref_NNs_x = basis_rotation(pref_NNs_x)
        PDFs_unnormalized = layer_evln(pref_NNs_x)

        return PDFs_unnormalized

    PDFs_unnormalized = compute_unnormalized_pdf(pdf_input)

    ## Normalization and sum rules, produces normalized PDF
    #   pdf_r(x) = updf_r(x) * Normalization(updf_r(integration_xgrid))
    # The normalization layer is shared across replicas (but evaluated at each replica separately)
    #
    if impose_sumrule:
        sumrule_layer, integrator_input = generate_msr_model_and_grid(
            fitbasis=fitbasis, mode=impose_sumrule, scaler=scaler, replica_seeds=all_seed
        )
        model_input["xgrid_integration"] = integrator_input

        # We need a second unnormalized PDF evaluated on the integrated grid
        PDFs_integration_grid = compute_unnormalized_pdf(integrator_input)

        # Photon contribution to the sum rule
        if photons:
            # add batch and flavor dimensions
            ph_tensor = op.numpy_to_tensor(photons.integral)
            photon_integrals = op.batchit(op.batchit(ph_tensor))
        else:
            photon_integrals = op.numpy_to_tensor(np.zeros((1, num_replicas, 1)))

        PDFs = sumrule_layer(
            {
                "pdf_x": PDFs_unnormalized,
                "pdf_xgrid_integration": PDFs_integration_grid,
                "xgrid_integration": integrator_input,
                # The photon is treated separately, need to get its integrals to normalize the pdf
                "photon_integral": photon_integrals,
            }
        )
    else:
        PDFs = PDFs_unnormalized
        sumrule_layer = lambda x: x

    ## Include the photon in the PDF for QED-enabled fits
    # (by default the entry corresponding to the photon is set to 0)
    if photons:
        layer_photon = AddPhoton(photons=photons, name="add_photon")
        PDFs = layer_photon(PDFs)

    # Return a PDF without a replica axis, to extract single replicas from an ensemble
    if not replica_axis:
        PDFs = Lambda(lambda pdfs: pdfs[:, 0], name="remove_replica_axis")(PDFs)

    return MetaModel(model_input, PDFs, name="PDFs", scaler=scaler)


# TODO: is there a way of keeping sincronized the input of this function and ReplicaSettings
# beyond a test of it? In principle we might want to have the arguments explicitly here...
def _generate_nn(
    input_layer: Input,
    replica_idx: int = 0,
    seed: int = None,
    nodes: list[int] = None,
    activations: list[str] = None,
    architecture: str = "dense",
    initializer: str = None,
    dropout_rate: float = 0.0,
    regularizer: str = None,
    regularizer_args: dict = field(default_factory=dict),
) -> MetaModel:
    """
    Create a Neural Network according to the input settings

    Parameters
    ----------
        input_layer: :py:class:`n3fit.backends.Input`
            input layer of the replica
        replica_idx: int
            Index of the replica used to name the PDF

        All other arguments follow exactly the documentation
        of ``ReplicaSettings``.
        See :py:class:`n3fit.model_gen.ReplicaSettings`


    Returns
    -------
        nn_pdf: MetaModel
            A single PDF NN model
    """
    reg = regularizer_selector(regularizer, **regularizer_args)
    *hidden_layers, n_flavours = nodes

    # Preparatory step: prepare a ``layer_generator`` function to iteratively create all layers
    # TODO: create a factory of layers instead of an ugly function
    # this layer generator takes the index of the layer (useful for seeding)
    # the output nodes of the layer
    # and the activation function

    if architecture == "dense_per_flavour":
        # Reset the last node in the list to be 1, we will then
        # repeat it n-times
        nodes = hidden_layers + [1]

        def layer_generator(i_layer, nodes_out, activation):
            """Generate the ``i_layer``-th dense_per_flavour layer for all replicas."""
            l_seed = int(seed + i_layer * n_flavours)
            initializers = [
                MetaLayer.select_initializer(initializer, seed=l_seed + b)
                for b in range(n_flavours)
            ]
            layer = base_layer_selector(
                architecture,
                kernel_initializer=initializers,
                units=int(nodes_out),
                activation=activation,
                basis_size=n_flavours,
            )
            return layer

    elif architecture == "dense":

        def layer_generator(i_layer, nodes_out, activation):
            kini = MetaLayer.select_initializer(initializer, seed=int(seed + i_layer))
            return base_layer_selector(
                architecture,
                kernel_initializer=kini,
                units=nodes_out,
                activation=activation,
                regularizer=reg,
            )

    else:
        raise ValueError(f"{architecture=} not recognized during model generation")

    # Use the previous layer generator to generate all layers
    previous_layer = input_layer
    for layer_idx, (nodes_out, activation) in enumerate(zip(nodes, activations)):
        layer = layer_generator(layer_idx, nodes_out, activation)

        # Apply the layer to the output of the previous one
        previous_layer = layer(previous_layer)

        # Add dropout if any to the second to last layer
        if dropout_rate > 0 and layer_idx == (len(hidden_layers) - 2):
            dropout_l = base_layer_selector("dropout", rate=dropout_rate)
            previous_layer = dropout_l(previous_layer)

    # In a dense-per-flavour, concatenate the last layer
    if architecture == "dense_per_flavour":
        concat = base_layer_selector("concatenate")
        previous_layer = concat(previous_layer)

    # Return the PDF model
    return MetaModel({"NN_input": input_layer}, previous_layer, name=f"{NN_PREFIX}_{replica_idx}")
