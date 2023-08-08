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

from n3fit.backends import Input, Lambda, MetaLayer, MetaModel, base_layer_selector
from n3fit.backends import operations as op
from n3fit.backends import regularizer_selector
from n3fit.layers import (
    DIS,
    DY,
    AddPhoton,
    FkRotation,
    FlavourToEvolution,
    ObsRotation,
    Preprocessing,
    losses,
)
from n3fit.layers.observable import is_unique
from n3fit.msr import generate_msr_model_and_grid
from validphys.photon.compute import Photon  # only used for type hint here


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
            loss = losses.LossInvcovmat(
                self.invcovmat, self.data, mask, covmat=self.covmat, name=self.name
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
                op_kwargs={"axis": 1},
                name=f"{self.name}_split",
            )
            sp_pdf = splitting_layer(pdf)
            output_layers = [obs(p) for obs, p in zip(self.observables, sp_pdf)]
        else:
            output_layers = [obs(pdf) for obs in self.observables]

        # Finally concatenate all observables (so that experiments are one single entitiy)
        ret = op.concatenate(output_layers, axis=2)
        if self.rotation is not None:
            ret = self.rotation(ret)
        return ret

    def __call__(self, pdf_layer, mask=None):
        loss_f = self._generate_loss(mask)
        experiment_prediction = self._generate_experimental_layer(pdf_layer)
        return loss_f(experiment_prediction)


def observable_generator(
    spec_dict, positivity_initial=1.0, integrability=False
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
        positivity_initial: float
            set the positivity lagrange multiplier for epoch 1

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
    model_obs_tr = []
    model_obs_vl = []
    model_obs_ex = []
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

        if spec_dict["positivity"]:
            # Positivity (and integrability, which is a special kind of positivity...)
            # enters only at the "training" part of the models
            obs_layer_tr = Obs_Layer(
                dataset.fktables_data,
                dataset.training_fktables(),
                operation_name,
                name=f"dat_{dataset_name}",
            )
            obs_layer_ex = obs_layer_vl = None
        elif spec_dict.get("data_transformation_tr") is not None:
            # Data transformation needs access to the full array of output data
            obs_layer_ex = Obs_Layer(
                dataset.fktables_data,
                dataset.fktables(),
                operation_name,
                name=f"exp_{dataset_name}",
            )
            obs_layer_tr = obs_layer_vl = obs_layer_ex
        else:
            obs_layer_tr = Obs_Layer(
                dataset.fktables_data,
                dataset.training_fktables(),
                operation_name,
                name=f"dat_{dataset_name}",
            )
            obs_layer_ex = Obs_Layer(
                dataset.fktables_data,
                dataset.fktables(),
                operation_name,
                name=f"exp_{dataset_name}",
            )
            obs_layer_vl = Obs_Layer(
                dataset.fktables_data,
                dataset.validation_fktables(),
                operation_name,
                name=f"val_{dataset_name}",
            )

        # If the observable layer found that all input grids are equal, the splitting will be None
        # otherwise the different xgrids need to be stored separately
        # Note: for pineappl grids, obs_layer_tr.splitting should always be None
        if obs_layer_tr.splitting is None:
            xgrid = dataset.fktables_data[0].xgrid
            model_inputs.append(xgrid)
            dataset_xsizes.append(len(xgrid))
        else:
            xgrids = [i.xgrid for i in dataset.fktables_data]
            model_inputs += xgrids
            dataset_xsizes.append(sum([len(i) for i in xgrids]))

        model_obs_tr.append(obs_layer_tr)
        model_obs_vl.append(obs_layer_vl)
        model_obs_ex.append(obs_layer_ex)

    # Check whether all xgrids of all observables in this experiment are equal
    # if so, simplify the model input
    if is_unique(model_inputs):
        model_inputs = model_inputs[0:1]
        dataset_xsizes = dataset_xsizes[0:1]

    # Reshape all inputs arrays to be (1, nx)
    model_inputs = np.concatenate(model_inputs).reshape(1, -1)

    full_nx = sum(dataset_xsizes)
    if spec_dict["positivity"]:
        out_positivity = ObservableWrapper(
            spec_name,
            model_obs_tr,
            dataset_xsizes,
            multiplier=positivity_initial,
            positivity=not integrability,
            integrability=integrability,
        )

        layer_info = {
            "inputs": model_inputs,
            "output_tr": out_positivity,
            "experiment_xsize": full_nx,
        }
        # For positivity we end here
        return layer_info

    # Generate the loss function and rotations of the final data (if any)
    if spec_dict.get("data_transformation_tr") is not None:
        obsrot_tr = ObsRotation(spec_dict.get("data_transformation_tr"))
        obsrot_vl = ObsRotation(spec_dict.get("data_transformation_vl"))
    else:
        obsrot_tr = None
        obsrot_vl = None

    out_tr = ObservableWrapper(
        spec_name,
        model_obs_tr,
        dataset_xsizes,
        invcovmat=spec_dict["invcovmat"],
        data=spec_dict["expdata"],
        rotation=obsrot_tr,
    )
    out_vl = ObservableWrapper(
        f"{spec_name}_val",
        model_obs_vl,
        dataset_xsizes,
        invcovmat=spec_dict["invcovmat_vl"],
        data=spec_dict["expdata_vl"],
        rotation=obsrot_vl,
    )
    out_exp = ObservableWrapper(
        f"{spec_name}_exp",
        model_obs_ex,
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
        "experiment_xsize": full_nx,
    }
    return layer_info


# Network generation functions
def generate_dense_network(
    nodes_in: int,
    nodes: int,
    activations: List[str],
    initializer_name: str = "glorot_normal",
    seed: int = 0,
    dropout_rate: float = 0.0,
    regularizer: str = None,
):
    """
    Generates a dense network

    the dropout rate, if selected, is set
    for the next to last layer (i.e., the last layer of the dense network before getting to
    the output layer for the basis choice)
    """
    list_of_pdf_layers = []
    number_of_layers = len(nodes)
    if dropout_rate > 0:
        dropout_layer = number_of_layers - 2
    else:
        dropout_layer = -1
    for i, (nodes_out, activation) in enumerate(zip(nodes, activations)):
        # if we have dropout set up, add it to the list
        if dropout_rate > 0 and i == dropout_layer:
            list_of_pdf_layers.append(base_layer_selector("dropout", rate=dropout_rate))

        # select the initializer and move the seed
        init = MetaLayer.select_initializer(initializer_name, seed=seed + i)

        # set the arguments that will define the layer
        arguments = {
            "kernel_initializer": init,
            "units": int(nodes_out),
            "activation": activation,
            "input_shape": (nodes_in,),
            "kernel_regularizer": regularizer,
        }

        layer = base_layer_selector("dense", **arguments)

        list_of_pdf_layers.append(layer)
        nodes_in = int(nodes_out)
    return list_of_pdf_layers


def generate_dense_per_flavour_network(
    nodes_in, nodes, activations, initializer_name="glorot_normal", seed=0, basis_size=8
):
    """
    For each flavour generates a dense network of the chosen size

    """
    list_of_pdf_layers = []
    number_of_layers = len(nodes)
    current_seed = seed
    for i, (nodes_out, activation) in enumerate(zip(nodes, activations)):
        initializers = []
        for _ in range(basis_size):
            # select the initializer and move the seed
            initializers.append(MetaLayer.select_initializer(initializer_name, seed=current_seed))
            current_seed += 1

        # set the arguments that will define the layer
        # but careful, the last layer must be nodes = 1
        # TODO the mismatch is due to the fact that basis_size
        # is set to the number of nodes of the last layer when it should
        # come from the runcard
        if i == number_of_layers - 1:
            nodes_out = 1
        arguments = {
            "kernel_initializer": initializers,
            "units": nodes_out,
            "activation": activation,
            "input_shape": (nodes_in,),
            "basis_size": basis_size,
        }

        layer = base_layer_selector("dense_per_flavour", **arguments)

        if i == number_of_layers - 1:
            # For the last layer, apply concatenate
            concat = base_layer_selector("concatenate")

            def output_layer(ilayer):
                result = layer(ilayer)
                return concat(result)

            list_of_pdf_layers.append(output_layer)
        else:
            list_of_pdf_layers.append(layer)

        nodes_in = int(nodes_out)
    return list_of_pdf_layers


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
    parallel_models: int = 1,
    photons: Photon = None,
    num_unique_As: int = 1,
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
    >>> pdf_model = pdfNN_layer_generator(nodes=[8], activations=['linear'], seed=[2,3], flav_info=fake_fl, parallel_models=2)

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
            whether to impose sumrules on the output pdf and which one to impose (All, MSR, VSR)
        scaler: callable
            Function to apply to the input. If given the input to the model
            will be a (1, None, 2) tensor where dim [:,:,0] is scaled
            When None, instead turn the x point into a (x, log(x)) pair
        parallel_models: int
            How many models should be trained in parallel
        photon: :py:class:`validphys.photon.compute.Photon`
            If given, gives the AddPhoton layer a function to compute a photon which will be added at the
            index 0 of the 14-size FK basis
            This same function will also be used to compute the MSR component for the photon
        num_unique_As: int
            Number of unique masses A, used for multi nuclei fits

    Returns
    -------
       pdf_models: list with a number equal to `parallel_models` of type n3fit.backends.MetaModel
            a model f(x) = y where x is a tensor (1, xgrid, 1) and y a tensor (1, xgrid, out)
    """
    # Parse the input configuration
    if seed is None:
        seed = parallel_models * [None]
    elif isinstance(seed, int):
        seed = parallel_models * [seed]

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

    number_of_layers = len(nodes)
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
        pdf_input_x = Input(shape=(None, pdf_input_dimensions), batch_size=1, name='scaledx_x')
        process_input = do_nothing
        extract_nn_input = Lambda(lambda x: op.op_gather_keep_dims(x, 0, axis=-1), name='x_scaled')
        extract_original = Lambda(lambda x: op.op_gather_keep_dims(x, 1, axis=-1), name='x')
    else:  # add log(x)
        pdf_input_x = Input(shape=(None, pdf_input_dimensions), batch_size=1, name='x')
        process_input = Lambda(lambda x: op.concatenate([x, op.op_log(x)], axis=-1), name='x_logx')
        extract_original = do_nothing
        extract_nn_input = do_nothing

    model_input = {"pdf_input_x": pdf_input_x}
    nn_input = {"NN_input_x": pdf_input_x}

    if num_unique_As > 1:
        pdf_input_A_stacked = Input(shape=(None,), batch_size=1, name='A_stacked')
        pdf_input_A_unique = Input(shape=(num_unique_As,), batch_size=1, name='A_unique')
        model_input["pdf_input_A_stacked"] = pdf_input_A_stacked
        nn_input["NN_input_A"] = pdf_input_A_stacked
        model_input["pdf_input_A_unique"] = pdf_input_A_unique

    if subtract_one:
        input_x_eq_1 = [1.0]
        if use_feature_scaling:
            input_x_eq_1 = scaler(input_x_eq_1)[0]
        # the layer that subtracts 1 from the NN output
        subtract_one_layer = Lambda(op.op_subtract, name='subtract_one')
        layer_x_eq_1 = op.numpy_to_input(np.array(input_x_eq_1).reshape(1, 1), name='x_eq_1')
        model_input["layer_x_eq_1"] = layer_x_eq_1

    # the layer that multiplies the NN output by the preprocessing factor
    apply_preprocessing_factor = Lambda(op.op_multiply, name='prefactor_times_NN')

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
            mode=impose_sumrule, scaler=scaler, photons=photons
        )
        model_input["integrator_input"] = integrator_input
        if num_unique_As == 1:
            nn_input_integration = {"NN_input_x": integrator_input}
        else:
            # Need to expand inputs to give each combination of integration gridpoint
            # and A value
            x_repeated, A_repeated = op.all_combinations(integrator_input, pdf_input_A_unique)
            nn_input_integration = {"NN_input_x": x_repeated, "NN_input_A": A_repeated}
    else:
        sumrule_layer = lambda x: x

    # Now we need a trainable network per replica to be trained in parallel
    pdf_models = []

    # Only these layers change from replica to replica:
    nn_replicas = []
    preprocessing_factor_replicas = []
    for i_replica, replica_seed in enumerate(seed):
        preprocessing_factor_replicas.append(
            Preprocessing(
                flav_info=flav_info,
                input_shape=(1,),
                name=f"preprocessing_factor_{i_replica}",
                seed=replica_seed + number_of_layers,
                large_x=not subtract_one,
            )
        )
        nn_replicas.append(
            generate_nn(
                layer_type=layer_type,
                input_dimensions=nn_input_dimensions,
                nodes=nodes,
                activations=activations,
                initializer_name=initializer_name,
                replica_seed=replica_seed,
                dropout=dropout,
                regularizer=regularizer,
                regularizer_args=regularizer_args,
                last_layer_nodes=last_layer_nodes,
                A_input=num_unique_As > 1,
                name=f"NN_{i_replica}",
            )
        )

    # All layers have been made, now we need to connect them,
    # do this in a function so we can call it for both grids and each replica
    # Since all layers are already made, they will be reused
    def compute_unnormalized_pdf(nn_input, neural_network, compute_preprocessing_factor):
        x = nn_input["NN_input_x"]
        # Preprocess the input grid
        x_nn_input = extract_nn_input(x)
        x_original = extract_original(x)
        x_processed = process_input(x_nn_input)

        # Compute the neural network output
        nn_input["NN_input_x"] = x_processed
        nn_output = neural_network(nn_input)
        if subtract_one:
            x_eq_1_processed = process_input(layer_x_eq_1)
            nn_input["NN_input_x"] = x_eq_1_processed
            nn_at_one = neural_network(nn_input)
            nn_output = subtract_one_layer([nn_output, nn_at_one])

        # Compute the preprocessing factor and multiply
        preprocessing_factor = compute_preprocessing_factor(x_original)
        pref_nn = apply_preprocessing_factor([nn_output, preprocessing_factor])

        # Apply basis rotation if needed
        if not basis_rotation.is_identity():
            pref_nn = basis_rotation(pref_nn)

        # Transform to FK basis
        pdf_unnormalized = layer_evln(pref_nn)

        return pdf_unnormalized

    # Finally compute the normalized PDFs for each replica
    pdf_models = []
    for i_replica, (preprocessing_factor, nn) in enumerate(
        zip(preprocessing_factor_replicas, nn_replicas)
    ):
        pdf_unnormalized = compute_unnormalized_pdf(nn_input, nn, preprocessing_factor)

        if impose_sumrule:
            pdf_integration_grid = compute_unnormalized_pdf(
                nn_input_integration, nn, preprocessing_factor
            )
            pdf_normalized = sumrule_layer(
                {
                    "pdf_x": pdf_unnormalized,
                    "pdf_xgrid_integration": pdf_integration_grid,
                    "xgrid_integration": integrator_input,
                    # The photon is treated separately, need to get its integrals to normalize the pdf
                    "photon_integral": op.numpy_to_tensor(
                        0.0 if not photons else photons.integral[i_replica]
                    ),
                }
            )
            pdf = pdf_normalized
        else:
            pdf = pdf_unnormalized

        if photons:
            # Add in the photon component
            pdf = layer_photon(pdf, i_replica)

        # Create the model
        pdf_model = MetaModel(model_input, pdf, name=f"PDF_{i_replica}", scaler=scaler)
        pdf_models.append(pdf_model)

    return pdf_models


def generate_nn(
    layer_type: str,
    input_dimensions: int,
    nodes: List[int],
    activations: List[str],
    initializer_name: str,
    replica_seed: int,
    dropout: float,
    regularizer: str,
    regularizer_args: dict,
    last_layer_nodes: int,
    name: str,
    A_input: bool,
) -> MetaModel:
    """
    Create the part of the model that contains all of the actual neural network
    layers.
    """
    common_args = {
        'nodes_in': input_dimensions,
        'nodes': nodes,
        'activations': activations,
        'initializer_name': initializer_name,
        'seed': replica_seed,
    }
    if layer_type == "dense":
        reg = regularizer_selector(regularizer, **regularizer_args)
        list_of_pdf_layers = generate_dense_network(
            **common_args, dropout_rate=dropout, regularizer=reg
        )
    elif layer_type == "dense_per_flavour":
        list_of_pdf_layers = generate_dense_per_flavour_network(
            **common_args, basis_size=last_layer_nodes
        )

    # Note: using a Sequential model would be more appropriate, but it would require
    # creating a MetaSequential model.
    x = Input(shape=(None, input_dimensions), batch_size=1, name='xgrids_processed')
    input_dict = {'NN_input_x': x}
    if A_input:
        A = Input(shape=(None, 1), batch_size=1, name='A')
        x = Concatenate(axis=-1)([x, A])
        input_dict['NN_input_A'] = A

    pdf = x
    for layer in list_of_pdf_layers:
        pdf = layer(pdf)

    model = MetaModel(input_dict, pdf, name=name)
    return model
