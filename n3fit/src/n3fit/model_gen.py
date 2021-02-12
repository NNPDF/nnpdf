"""
    Library of functions which generate the NN objects

    Contains:
        # observable_generator:
            Generates the output layers as functions
        # pdfNN_layer_generator:
            Generates the PDF NN layer to be fitted
"""
import numpy as np
import n3fit.msr as msr_constraints
from n3fit.layers import DIS, DY, Mask, ObsRotation
from n3fit.layers import Preprocessing, FkRotation, FlavourToEvolution

from n3fit.backends import MetaModel, Input
from n3fit.backends import operations
from n3fit.backends import losses
from n3fit.backends import MetaLayer, Lambda
from n3fit.backends import base_layer_selector, regularizer_selector


def observable_generator(spec_dict, positivity_initial=1.0, integrability=False):  # pylint: disable=too-many-locals
    """
    This function generates the observable model for each experiment.
    These are models which takes as input a PDF tensor (1 x size_of_xgrid x flavours) and outputs
    the result of the observable for each contained dataset (n_points,)

    An experiment contains an fktable, which is loaded by the convolution layer
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
            - `loss` : loss function (unmasked)
            - `output_tr`: output layer (training)
            - `loss_tr` : loss function (training)
            - `output_vl`: output layer (validation)
            - `loss_vl` : loss function (validation)
    """
    spec_name = spec_dict["name"]
    dataset_xsizes = []
    model_obs_tr = []
    model_obs_vl = []
    model_obs_ex = []
    model_inputs = []
    # The first step is to compute the observable for each of the datasets
    for dataset_dict in spec_dict["datasets"]:
        # Get the generic information of the dataset
        dataset_name = dataset_dict["name"]

        # Look at what kind of layer do we need for this dataset
        if dataset_dict["hadronic"]:
            Obs_Layer = DY
        else:
            Obs_Layer = DIS

        # Set the operation (if any) to be applied to the fktables of this dataset
        operation_name = dataset_dict["operation"]

        # Now generate the observable layer, which takes the following information:
        # operation name
        # dataset name
        # list of fktable_dictionaries
        #   these will then be used to check how many different pdf inputs are needed
        #   (and convolutions if given the case)
        obs_layer_tr = Obs_Layer(
            dataset_dict["fktables"],
            dataset_dict["tr_fktables"],
            operation_name,
            name=f"dat_{dataset_name}",
        )
        if spec_dict["positivity"]:
            obs_layer_ex = obs_layer_tr
            obs_layer_vl = obs_layer_tr
        else:
            obs_layer_ex = Obs_Layer(
                dataset_dict["fktables"],
                dataset_dict["ex_fktables"],
                operation_name,
                name=f"exp_{dataset_name}",
            )
            obs_layer_vl = Obs_Layer(
                dataset_dict["fktables"],
                dataset_dict["vl_fktables"],
                operation_name,
                name=f"val_{dataset_name}",
            )

        # Data transformation might need access to the full array of output data
        # therefore the validation and training layers should point to the full exp
        if spec_dict.get("data_transformation") is not None:
            obs_layer_tr = obs_layer_ex
            obs_layer_vl = obs_layer_ex

        # To know how many xpoints we compute we are duplicating functionality from obs_layer
        if obs_layer_tr.splitting is None:
            xgrid = dataset_dict["fktables"][0]["xgrid"]
            model_inputs.append(xgrid)
            dataset_xsizes.append(xgrid.shape[1])
        else:
            xgrids = [i["xgrid"] for i in dataset_dict["fktables"]]
            model_inputs += xgrids
            dataset_xsizes.append(sum([i.shape[1] for i in xgrids]))

        model_obs_tr.append(obs_layer_tr)
        model_obs_vl.append(obs_layer_vl)
        model_obs_ex.append(obs_layer_ex)

    # Prepare a concatenation as experiments are one single entity formed by many datasets
    def gen_concat(name):
        return operations.as_layer(operations.concatenate, op_kwargs={"axis": 1}, name=name)

    # Tensorflow operations have ugly name,
    # we want the final observables to be named just {spec_name} (with'val/exp' if needed)
    tr_name = spec_name
    vl_name = f"{spec_name}_val"
    ex_name = f"{spec_name}_exp"
    concat_ex = gen_concat(ex_name)
    # For data transformation all concatenations are the same
    if spec_dict.get("data_transformation") is None:
        concat_tr = gen_concat(tr_name)
        concat_vl = gen_concat(vl_name)
    else:
        concat_tr = concat_ex
        concat_vl = concat_ex

    # creating the experiment as a model turns out to bad for performance
    def experiment_layer(pdf, model_obs=model_obs_ex, concat=concat_ex, rotation=None, datasets_out=None):
        """ By default works with the experiment observable """
        output_layers = []
        # First split the pdf layer into the different datasets if needed
        if len(dataset_xsizes) > 1:
            splitting_layer = operations.as_layer(
                operations.split,
                op_args=[dataset_xsizes],
                op_kwargs={"axis": 1},
                name=f"{spec_name}_split",
            )
            split_pdf = splitting_layer(pdf)
        else:
            split_pdf = [pdf]
        # every obs gets its share of the split
        for partial_pdf, obs in zip(split_pdf, model_obs):
            obs_output = obs(partial_pdf)
            if datasets_out and obs.name[4:] in datasets_out:
                mask_out = Mask(c=0.0, name=f"zero_{obs.name}")
                obs_output = mask_out(obs_output)
            output_layers.append(obs_output)
        # Concatenate all datasets as experiments are one single entity if needed
        ret = concat(output_layers)
        if rotation is not None:
            ret = rotation(ret)
        return ret

    # Now create the model for this experiment
    full_nx = sum(dataset_xsizes)

    if spec_dict["positivity"]:
        out_mask = Mask(
            c=positivity_initial,
            axis=1,
            name=spec_name,
        )

        def out_positivity(pdf_layer, datasets_out=None):
            exp_result = experiment_layer(pdf_layer)
            return out_mask(exp_result)

        if integrability:
            loss = losses.l_integrability()
        else:
            loss = losses.l_positivity()

        layer_info = {
            "inputs": model_inputs,
            "output_tr": out_positivity,
            "loss_tr": loss,
            "experiment_xsize": full_nx,
        }
        return layer_info

    invcovmat_tr = spec_dict["invcovmat"]
    invcovmat_vl = spec_dict["invcovmat_vl"]
    invcovmat = spec_dict["invcovmat_true"]

    # Generate the loss function and rotations of the final data (if any)
    if spec_dict.get("data_transformation") is not None:
        # The rotation is the last layer so it should carry The Name
        obsrot_tr = ObsRotation(spec_dict.get("data_transformation"), name=tr_name)
        obsrot_vl = ObsRotation(spec_dict.get("data_transformation_vl"), name=vl_name)
        loss_tr = losses.l_diaginvcovmat(invcovmat_tr)
        loss_vl = losses.l_diaginvcovmat(invcovmat_vl)
    else:
        obsrot_tr = None
        obsrot_vl = None
        loss_tr = losses.l_invcovmat(invcovmat_tr)
        # TODO At this point we need to intercept the data and compile the loss with it
        # then the validation must have a list of None as an output
        loss_vl = losses.l_invcovmat(invcovmat_vl)
    loss = losses.l_invcovmat(invcovmat)

    def out_tr(pdf_layer, datasets_out=None):
        exp_result = experiment_layer(
            pdf_layer, model_obs=model_obs_tr, concat=concat_tr, datasets_out=datasets_out, rotation=obsrot_tr
        )
        return exp_result

    def out_vl(pdf_layer, datasets_out=None):
        exp_result = experiment_layer(
            pdf_layer, model_obs=model_obs_vl, concat=concat_vl, datasets_out=datasets_out, rotation=obsrot_vl
        )
        return exp_result

    layer_info = {
        "inputs": model_inputs,
        "output": experiment_layer,
        "loss": loss,
        "output_tr": out_tr,
        "loss_tr": loss_tr,
        "output_vl": out_vl,
        "loss_vl": loss_vl,
        "experiment_xsize": full_nx,
    }

    return layer_info


# Network generation functions
def generate_dense_network(
    nodes_in,
    nodes,
    activations,
    initializer_name="glorot_normal",
    seed=0,
    dropout_rate=0.0,
    regularizer=None,
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
    inp=2,
    nodes=None,
    activations=None,
    initializer_name="glorot_normal",
    layer_type="dense",
    flav_info=None,
    fitbasis="NN31IC",
    out=14,
    seed=None,
    dropout=0.0,
    regularizer=None,
    regularizer_args=None,
    impose_sumrule=False,
    scaler=None,
):  # pylint: disable=too-many-locals
    """
    Generates the PDF model which takes as input a point in x (from 0 to 1)
    and outputs a basis of 14 PDFs.
    It generates the preprocessing of the x into a set (x, log(x)),
    the arbitrary NN to fit the form of the PDF
    and the preprocessing factors.

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

    4. Create a preprocessing layer (that takes as input the same tensor x as the NN)
    and multiply it to the NN. We have now:
        N(x)_{j} * x^{1-alpha_{j}} * (1-x)^{beta_{j}}

    5. Create a rotation layer (that takes as input the previous step) and applies it
    so that we get f(x)_{i}

    Finally we output the final answer as well as the list of all generating functions
    in the model for easy usage within `n3fit`.

    Parameters
    ----------
        inp: int
            dimension of the xgrid. If inp=2, turns the x point into a (x, log(x)) pair
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
        seed: int
            seed to initialize the NN
        dropout: float
            rate of dropout layer by layer
        impose_sumrule: bool
            whether to impose sumrule on the output pdf model

    Returns
    -------
        model_pdf: n3fit.backends.MetaModel
            a model f(x) = y where x is a tensor (1, xgrid, 1) and y a tensor (1, xgrid, out)
    """
    if nodes is None:
        nodes = [15, 8]
    ln = len(nodes)

    if scaler:
        inp = 1

    if activations is None:
        activations = ["tanh", "linear"]
    elif callable(activations):
        # hyperopt passes down a function to generate dynamically the list of
        # activations functions
        activations = activations(ln)

    if regularizer_args is None:
        regularizer_args = dict()
    # Safety check
    number_of_layers = len(nodes)
    number_of_activations = len(activations)
    if number_of_layers != number_of_activations:
        raise ValueError(
            "Number of activation functions does not match number of layers @ model_gen.py"
        )
    # The number of nodes in the last layer is equal to the number of fitted flavours (== len(flav_info))
    last_layer_nodes = nodes[-1]

    if layer_type == "dense":
        reg = regularizer_selector(regularizer, **regularizer_args)
        list_of_pdf_layers = generate_dense_network(
            inp,
            nodes,
            activations,
            initializer_name,
            seed=seed,
            dropout_rate=dropout,
            regularizer=reg,
        )
    elif layer_type == "dense_per_flavour":
        # Define the basis size attending to the last layer in the network
        # TODO: this information should come from the basis information
        #       once the basis information is passed to this class
        list_of_pdf_layers = generate_dense_per_flavour_network(
            inp, nodes, activations, initializer_name, seed=seed, basis_size=last_layer_nodes,
        )

    placeholder_input = Input(shape=(None, 1), batch_size=1)
    subtract_one = False
    process_input = Lambda(lambda x: x)
    input_x_eq_1 = [1.0]
    # When scaler is active we also want to do the subtraction of large x
    # TODO: make it its own option (i.e., one could want to use this without using scaler)
    if scaler:
        # change the input domain [0,1] -> [-1,1]
        process_input = Lambda(lambda  x: 2*x-1)
        subtract_one = True
        input_x_eq_1 = scaler([1.0])[0]
        placeholder_input = Input(shape=(None, 2), batch_size=1)
    elif inp==2:
        # If the input is of type (x, logx)
        # create a x --> (x, logx) layer to preppend to everything
        process_input = Lambda(lambda x: operations.concatenate([x, operations.op_log(x)], axis=-1))

    model_input = [placeholder_input]
    if subtract_one:
        layer_x_eq_1 = operations.numpy_to_input(np.array(input_x_eq_1).reshape(1,1))
        model_input.append(layer_x_eq_1)

    def dense_me(x):
        """Takes an input tensor `x` and applies all layers
        from the `list_of_pdf_layers` in order"""
        processed_x = process_input(x)
        curr_fun = list_of_pdf_layers[0](processed_x)

        for dense_layer in list_of_pdf_layers[1:]:
            curr_fun = dense_layer(curr_fun)
        return curr_fun

    # Preprocessing layer (will be multiplied to the last of the denses)
    preproseed = seed + number_of_layers
    layer_preproc = Preprocessing(
        input_shape=(1,),
        name="pdf_prepro",
        flav_info=flav_info,
        seed=preproseed,
        output_dim=last_layer_nodes,
        large_x = not subtract_one
    )
    # Basis rotation
    basis_rotation = FlavourToEvolution(flav_info=flav_info, fitbasis=fitbasis)

    # Evolution layer
    layer_evln = FkRotation(input_shape=(last_layer_nodes,), output_dim=out)


    # Apply extrapolation and basis
    def layer_fitbasis(x):
        """ The tensor x has a expected shape of (1, None, {1,2})
        where x[...,0] corresponds to the feature_scaled input and x[...,-1] the original input
        """
        x_scaled = operations.op_gather_keep_dims(x, 0, axis=-1)
        x_original = operations.op_gather_keep_dims(x, -1, axis=-1)

        nn_output = dense_me(x_scaled)
        if subtract_one:
            nn_at_one = dense_me(layer_x_eq_1)
            nn_output = operations.op_subtract([nn_output, nn_at_one])

        ret = operations.op_multiply([nn_output, layer_preproc(x_original)])
        if basis_rotation.is_identity():
            # if we don't need to rotate basis we don't want spurious layers
            return ret
        return basis_rotation(ret)

    # Rotation layer, changes from the 8-basis to the 14-basis
    def layer_pdf(x):
        return layer_evln(layer_fitbasis(x))

    # Impose sumrule if necessary
    if impose_sumrule:
        layer_pdf, integrator_input = msr_constraints.msr_impose(layer_fitbasis, layer_pdf, scaler=scaler)
        model_input.append(integrator_input)
    else:
        integrator_input = None

    pdf_model = MetaModel(model_input, layer_pdf(placeholder_input), name="PDF", scaler=scaler)

    return pdf_model
