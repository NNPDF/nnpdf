"""
    Library of functions which generate the NN objects

    Contains:
        # observable_generator:
            Generates the output layers as functions
        # pdfNN_layer_generator:
            Generates the PDF NN layer to be fitted


"""

from dataclasses import dataclass
from typing import Callable, List

import numpy as np

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
from n3fit.backends import regularizer_selector
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
    invcovmat: np.array = None
    covmat: np.array = None
    multiplier: float = 1.0
    integrability: bool = False
    positivity: bool = False
    data: np.array = None
    rotation: ObsRotation = None  # only used for diagonal covmat

    def _generate_loss(self, mask=None):
        """Generates the corresponding loss function depending on the values the wrapper
        was initialized with"""
        if self.invcovmat is not None:
            if self.rotation:
                # If we have a matrix diagonal only, padd with 0s and hope it's not too heavy on memory
                invcovmat_matrix = (
                    np.eye(self.invcovmat.shape[-1]) * self.invcovmat[..., np.newaxis]
                )
                if self.covmat is not None:
                    covmat_matrix = np.eye(self.covmat.shape[-1]) * self.covmat[..., np.newaxis]
                else:
                    covmat_matrix = self.covmat
            else:
                covmat_matrix = self.covmat
                invcovmat_matrix = self.invcovmat
            loss = losses.LossInvcovmat(
                invcovmat_matrix, self.data, mask, covmat=covmat_matrix, name=self.name
            )
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
            splitting_layer = op.as_layer(
                op.split,
                op_args=[self.dataset_xsizes],
                op_kwargs={"axis": 2},
                name=f"{self.name}_split",
            )
            sp_pdf = splitting_layer(pdf)
            output_layers = [obs(p) for obs, p in zip(self.observables, sp_pdf)]
        else:
            output_layers = [obs(pdf) for obs in self.observables]

        # Finally concatenate all observables (so that experiments are one single entity)
        ret = op.concatenate(output_layers, axis=-1)

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
    mask_array=None,
    training_data=None,
    validation_data=None,
    invcovmat_tr=None,
    invcovmat_vl=None,
    positivity_initial=1.0,
    integrability=False,
    n_replicas=1,
):  # pylint: disable=too-many-locals
    """
    This function generates the observable models for each experiment.
    These are models which takes as input a PDF tensor (1 x size_of_xgrid x flavours) and outputs
    the result of the observable for each contained dataset (n_points,).

    In summary the model has the following structure:
        One experiment layer, made of any number of observable layers.
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
    if mask_array is not None:
        tr_mask_layer = Mask(mask_array, name=f"trmask_{spec_name}")
        vl_mask_layer = Mask(~mask_array, name=f"vlmask_{spec_name}")
    else:
        tr_mask_layer = None
        vl_mask_layer = None

    # Make rotations of the final data (if any)
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
        invcovmat=invcovmat_tr,
        data=training_data,
        rotation=obsrot,
    )
    out_vl = ObservableWrapper(
        f"{spec_name}_val",
        model_observables,
        vl_mask_layer,
        dataset_xsizes,
        invcovmat=invcovmat_vl,
        data=validation_data,
        rotation=obsrot,
    )
    out_exp = ObservableWrapper(
        f"{spec_name}_exp",
        model_observables,
        None,
        dataset_xsizes,
        invcovmat=spec_dict["invcovmat_true"],
        covmat=spec_dict["covmat"],
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


def generate_pdf_model(
    nodes: List[int] = None,
    activations: List[str] = None,
    initializer_name: str = "glorot_normal",
    layer_type: str = "dense",
    flav_info: dict = None,
    fitbasis: str = "NN31IC",
    out: int = 14,
    seed: int = None,
    dropout: float = 0.0,
    regularizer: str = None,
    regularizer_args: dict = None,
    impose_sumrule: str = None,
    scaler: Callable = None,
    num_replicas: int = 1,
    photons: Photon = None,
):
    """
    Wrapper around pdfNN_layer_generator to allow the generation of single replica models.

    Parameters:
    -----------
        see model_gen.pdfNN_layer_generator

    Returns
    -------
        pdf_model: MetaModel
            pdf model, with `single_replica_generator` attached in a list as an attribute
    """
    joint_args = {
        "nodes": nodes,
        "activations": activations,
        "initializer_name": initializer_name,
        "layer_type": layer_type,
        "flav_info": flav_info,
        "fitbasis": fitbasis,
        "out": out,
        "dropout": dropout,
        "regularizer": regularizer,
        "regularizer_args": regularizer_args,
        "impose_sumrule": impose_sumrule,
        "scaler": scaler,
    }

    pdf_model = pdfNN_layer_generator(
        **joint_args, seed=seed, num_replicas=num_replicas, photons=photons
    )

    # Note that the photons are passed unchanged to the single replica generator
    # computing the photon requires running fiatlux which takes 30' per replica
    # and so at the moment parallel photons are disabled with a check in checks.py
    # In order to enable it `single_replica_generator` must take the index of the replica
    # to select the appropiate photon as all of them will be computed and fixed before the fit

    # this is necessary to be able to convert back to single replica models after training
    single_replica_generator = lambda: pdfNN_layer_generator(
        **joint_args, seed=0, num_replicas=1, photons=photons, replica_axis=False
    )
    pdf_model.single_replica_generator = single_replica_generator

    return pdf_model


def pdfNN_layer_generator(
    nodes: List[int] = None,
    activations: List[str] = None,
    initializer_name: str = "glorot_normal",
    layer_type: str = "dense",
    flav_info: dict = None,
    fitbasis: str = "NN31IC",
    out: int = 14,
    seed: int = None,
    dropout: float = 0.0,
    regularizer: str = None,
    regularizer_args: dict = None,
    impose_sumrule: str = None,
    scaler: Callable = None,
    num_replicas: int = 1,
    photons: Photon = None,
    replica_axis: bool = True,
):  # pylint: disable=too-many-locals
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
    >>> from n3fit.model_gen import pdfNN_layer_generator
    >>> from validphys.pdfgrids import xplotting_grid
    >>> fake_fl = [{'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'cbar', 's', 'sbar']]
    >>> fake_x = np.linspace(1e-3,0.8,3)
    >>> pdf_model = pdfNN_layer_generator(nodes=[8], activations=['linear'], seed=[2,3], flav_info=fake_fl, num_replicas=2)

    Parameters
    ----------
        nodes: list(int)
            list of the number of nodes per layer of the PDF NN. Default: [15,8]
        activation: list
            list of activation functions to apply to each layer. Default: ["tanh", "linear"]
            if the number of activation function does not match the number of layers, it will add
            copies of the first activation function found
        initializer_name: str
            selects the initializer of the weights of the NN. Default: glorot_normal
        layer_type: str
            selects the type of architecture of the NN. Default: dense
        flav_info: dict
            dictionary containing the information about each PDF (basis dictionary in the runcard)
            to be used by Preprocessing
        out: int
            number of output flavours of the model (default 14)
        seed: list(int)
            seed to initialize the NN
        dropout: float
            rate of dropout layer by layer
        impose_sumrule: str
            whether to impose sumrules on the output pdf and which one to impose (All, MSR, VSR, TSR)
        scaler: callable
            Function to apply to the input. If given the input to the model
            will be a (1, None, 2) tensor where dim [:,:,0] is scaled
            When None, instead turn the x point into a (x, log(x)) pair
        num_replicas: int
            How many models should be trained in parallel
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
    # Parse the input configuration
    if seed is None:
        seed = num_replicas * [None]
    elif isinstance(seed, int):
        seed = num_replicas * [seed]

    if nodes is None:
        nodes = [15, 8]
    ln = len(nodes)

    if impose_sumrule is None:
        impose_sumrule = "All"

    if activations is None:
        activations = ["tanh", "linear"]
    elif callable(activations):
        # hyperopt passes down a function to generate dynamically the list of
        # activations functions
        activations = activations(ln)

    if regularizer_args is None:
        regularizer_args = dict()

    # The number of nodes in the last layer is equal to the number of fitted flavours
    last_layer_nodes = nodes[-1]  # (== len(flav_info))

    # Process input options. There are 2 options:
    # 1. Scale the input
    # 2. Concatenate log(x) to the input
    use_feature_scaling = scaler is not None

    # When scaler is active we also want to do the subtraction of large x
    # TODO: make it its own option (i.e., one could want to use this without using scaler)
    subtract_one = use_feature_scaling

    # Feature scaling happens before the pdf model and changes x->(scaler(x), x),
    # so it adds an input dimension
    pdf_input_dimensions = 2 if use_feature_scaling else 1
    # Adding of logs happens inside, but before the NN and adds a dimension there
    nn_input_dimensions = 1 if use_feature_scaling else 2

    # Define the main input
    do_nothing = lambda x: x
    if use_feature_scaling:
        pdf_input = Input(shape=(None, pdf_input_dimensions), batch_size=1, name="scaledx_x")
        process_input = do_nothing
        extract_nn_input = Lambda(lambda x: op.op_gather_keep_dims(x, 0, axis=-1), name="x_scaled")
        extract_original = Lambda(lambda x: op.op_gather_keep_dims(x, 1, axis=-1), name="pdf_input")
    else:  # add log(x)
        pdf_input = Input(shape=(None, pdf_input_dimensions), batch_size=1, name="pdf_input")
        process_input = Lambda(lambda x: op.concatenate([x, op.op_log(x)], axis=-1), name="x_logx")
        extract_original = do_nothing
        extract_nn_input = do_nothing

    model_input = {"pdf_input": pdf_input}

    if subtract_one:
        input_x_eq_1 = [1.0]
        if use_feature_scaling:
            input_x_eq_1 = scaler(input_x_eq_1)[0]
        # the layer that subtracts 1 from the NN output
        subtract_one_layer = Lambda(op.op_subtract, name="subtract_one")
        layer_x_eq_1 = op.numpy_to_input(np.array(input_x_eq_1).reshape(1, 1), name="x_eq_1")
        model_input["layer_x_eq_1"] = layer_x_eq_1

    # the layer that multiplies the NN output by the preprocessing factor
    apply_preprocessing_factor = Lambda(op.op_multiply, name="prefactor_times_NN")

    # Photon layer
    layer_photon = AddPhoton(photons=photons, name="add_photon")

    # Basis rotation
    basis_rotation = FlavourToEvolution(
        flav_info=flav_info, fitbasis=fitbasis, name="pdf_evolution_basis"
    )

    # Evolution layer
    layer_evln = FkRotation(input_shape=(last_layer_nodes,), output_dim=out, name="pdf_FK_basis")

    # Normalization and sum rules
    if impose_sumrule:
        sumrule_layer, integrator_input = generate_msr_model_and_grid(
            fitbasis=fitbasis, mode=impose_sumrule, scaler=scaler, replica_seeds=seed
        )
        model_input["xgrid_integration"] = integrator_input
    else:
        sumrule_layer = lambda x: x

    compute_preprocessing_factor = Preprocessing(
        flav_info=flav_info,
        input_shape=(1,),
        name=PREPROCESSING_LAYER_ALL_REPLICAS,
        replica_seeds=seed,
        large_x=not subtract_one,
    )

    nn_replicas = generate_nn(
        layer_type=layer_type,
        nodes_in=nn_input_dimensions,
        nodes=nodes,
        activations=activations,
        initializer_name=initializer_name,
        replica_seeds=seed,
        dropout=dropout,
        regularizer=regularizer,
        regularizer_args=regularizer_args,
        last_layer_nodes=last_layer_nodes,
    )

    # The NN subtracted by NN(1), if applicable
    def nn_subtracted(x):
        NNs_x = nn_replicas(x)

        if subtract_one:
            x_eq_1_processed = process_input(layer_x_eq_1)
            NNs_x_1 = nn_replicas(x_eq_1_processed)
            NNs_x = subtract_one_layer([NNs_x, NNs_x_1])

        return NNs_x

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

        # Apply basis rotation if needed
        if not basis_rotation.is_identity():
            pref_NNs_x = basis_rotation(pref_NNs_x)

        # Transform to FK basis
        PDFs_unnormalized = layer_evln(pref_NNs_x)

        return PDFs_unnormalized

    PDFs_unnormalized = compute_unnormalized_pdf(pdf_input)

    if impose_sumrule:
        PDFs_integration_grid = compute_unnormalized_pdf(integrator_input)

        if photons:
            # add batch and flavor dimensions
            ph_tensor = op.numpy_to_tensor(photons.integral)
            photon_integrals = op.batchit(op.batchit(ph_tensor))
        else:
            photon_integrals = op.numpy_to_tensor(np.zeros((1, num_replicas, 1)))

        PDFs_normalized = sumrule_layer(
            {
                "pdf_x": PDFs_unnormalized,
                "pdf_xgrid_integration": PDFs_integration_grid,
                "xgrid_integration": integrator_input,
                # The photon is treated separately, need to get its integrals to normalize the pdf
                "photon_integral": photon_integrals,
            }
        )
        PDFs = PDFs_normalized
    else:
        PDFs = PDFs_unnormalized

    if photons:
        PDFs = layer_photon(PDFs)

    if not replica_axis:
        PDFs = Lambda(lambda pdfs: pdfs[:, 0], name="remove_replica_axis")(PDFs)

    pdf_model = MetaModel(model_input, PDFs, name=f"PDFs", scaler=scaler)
    return pdf_model


def generate_nn(
    layer_type: str,
    nodes_in: int,
    nodes: List[int],
    activations: List[str],
    initializer_name: str,
    replica_seeds: List[int],
    dropout: float,
    regularizer: str,
    regularizer_args: dict,
    last_layer_nodes: int,
) -> MetaModel:
    """
    Create the part of the model that contains all of the actual neural network
    layers, for each replica.

    Parameters
    ----------
        layer_type: str
            Type of layer to use. Can be "dense" or "dense_per_flavour".
        nodes_in: int
            Number of nodes in the input layer.
        nodes: List[int]
            Number of nodes in each hidden layer.
        activations: List[str]
            Activation function to use in each hidden layer.
        initializer_name: str
            Name of the initializer to use.
        replica_seeds: List[int]
            List of seeds to use for each replica.
        dropout: float
            Dropout rate to use (if 0, no dropout is used).
        regularizer: str
            Name of the regularizer to use.
        regularizer_args: dict
            Arguments to pass to the regularizer.
        last_layer_nodes: int
            Number of nodes in the last layer.

    Returns
    -------
        nn_replicas: MetaModel
            Single model containing all replicas.
    """
    nodes_list = list(nodes)  # so we can modify it
    x_input = Input(shape=(None, nodes_in), batch_size=1, name="NN_input")
    reg = regularizer_selector(regularizer, **regularizer_args)

    if layer_type == "dense_per_flavour":
        # set the arguments that will define the layer
        # but careful, the last layer must be nodes = 1
        # TODO the mismatch is due to the fact that basis_size
        # is set to the number of nodes of the last layer when it should
        # come from the runcard
        nodes_list[-1] = 1
        basis_size = last_layer_nodes

        def layer_generator(i_layer, nodes_out, activation):
            """Generate the ``i_layer``-th dense_per_flavour layer for all replicas."""
            layers = []
            for replica_seed in replica_seeds:
                seed = replica_seed + i_layer * basis_size
                initializers = [
                    MetaLayer.select_initializer(initializer_name, seed=seed + b)
                    for b in range(basis_size)
                ]
                layer = base_layer_selector(
                    layer_type,
                    kernel_initializer=initializers,
                    units=int(nodes_out),
                    activation=activation,
                    input_shape=(nodes_in,),
                    basis_size=basis_size,
                )
                layers.append(layer)

            return layers

    elif layer_type == "single_dense":

        # The checks should've triggered, but better safe than sorry
        if len(replica_seeds) > 1:
            raise ValueError("`single_dense` only valid with one replica")
        seed = replica_seeds[0]

        def layer_generator(i_layer, nodes_out, activation):
            return base_layer_selector(
                layer_type,
                kernel_initializer=MetaLayer.select_initializer(
                    initializer_name, seed=seed + i_layer
                ),
                units=nodes_out,
                activation=activation,
                input_shape=(nodes_in,),
                regularizer=reg,
            )

    elif layer_type == "dense":

        def layer_generator(i_layer, nodes_out, activation):
            """Generate the ``i_layer``-th MetaLayer.MultiDense layer for all replicas."""
            return base_layer_selector(
                layer_type,
                replica_seeds=replica_seeds,
                kernel_initializer=MetaLayer.select_initializer(initializer_name),
                base_seed=i_layer,
                units=int(nodes_out),
                activation=activation,
                is_first_layer=(i_layer == 0),
                regularizer=reg,
                name=f"multi_dense_{i_layer}",
            )

    else:
        raise ValueError(f"{layer_type=} not recognized during model generation")

    # First create all the layers
    # list_of_pdf_layers[d][r] is the layer at depth d for replica r
    list_of_pdf_layers = []
    for i_layer, (nodes_out, activation) in enumerate(zip(nodes_list, activations)):
        layers = layer_generator(i_layer, nodes_out, activation)
        list_of_pdf_layers.append(layers)
        nodes_in = int(nodes_out)

    # add dropout as second to last layer
    if dropout > 0:
        dropout_layer = base_layer_selector("dropout", rate=dropout)
        list_of_pdf_layers.insert(-2, dropout_layer)

    # In case of per flavour network, concatenate at the last layer
    if layer_type == "dense_per_flavour":
        concat = base_layer_selector("concatenate")
        list_of_pdf_layers[-1] = [lambda x: concat(layer(x)) for layer in list_of_pdf_layers[-1]]

    # Apply all layers to the input to create the models
    if layer_type in ("dense", "single_dense"):
        pdfs = x_input
        for layer in list_of_pdf_layers:
            pdfs = layer(pdfs)
        model = MetaModel({'NN_input': x_input}, pdfs, name=NN_LAYER_ALL_REPLICAS)

        return model

    pdfs = [layer(x_input) for layer in list_of_pdf_layers[0]]

    for layers in list_of_pdf_layers[1:]:
        # Since some layers (dropout) are shared, we have to treat them separately
        if type(layers) is list:
            pdfs = [layer(x) for layer, x in zip(layers, pdfs)]
        else:
            pdfs = [layers(x) for x in pdfs]

    # Wrap the pdfs in a MetaModel to enable getting/setting of weights later
    pdfs = [
        MetaModel({'NN_input': x_input}, pdf, name=f"{NN_PREFIX}_{i_replica}")(x_input)
        for i_replica, pdf in enumerate(pdfs)
    ]
    pdfs = Lambda(lambda nns: op.stack(nns, axis=1), name=f"stack_replicas")(pdfs)
    model = MetaModel({'NN_input': x_input}, pdfs, name=NN_LAYER_ALL_REPLICAS)

    return model
