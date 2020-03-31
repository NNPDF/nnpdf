"""
    Library of functions which generate the NN objects

    Contains:
        # observable_generator:
            Generates the output layers
        # pdfNN_layer_generator:
            Generates the PDF NN layer to be fitted
"""
from n3fit.layers import DIS
from n3fit.layers import DY
from n3fit.layers import Mask
from n3fit.layers import Preprocessing, Rotation

from n3fit.backends import operations
from n3fit.backends import losses
from n3fit.backends import MetaLayer
from n3fit.backends import base_layer_selector, concatenate, Lambda


def observable_generator(
    spec_dict, positivity_initial=None, positivity_multiplier=1.05, positivity_steps=300
):  # pylint: disable=too-many-locals
    """
    This function generates the observable generator.
    It loads the fktable in the convolution layer (hadronic or DIS) and prepares the loss function.
    If the dataset is a positivity dataset acts in consequence.

    Parameters
    ----------
        `spec_dict`
            a dictionary-like object containing the information of the observable
        `positivity_initial`
            if given, set this number as the positivity multiplier for epoch 1
        `positivity_multiplier`
            how much the positivity increases every 100 steps
        `positivity_steps`
            if positivity_initial is not given, computes the initial by assuming we want,
            after 100**positivity_steps epochs, to have the lambda of the runcard

    Returns
    ------
        'layer_info` a dictionary with:
            - `inputs`: input layer
            - `output`: output layer (unmasked)
            - `loss` : loss function (unmasked)
            - `output_tr`: output layer (training)
            - `loss_tr` : loss function (training)
            - `output_vl`: output layer (validation)
            - `loss_vl` : loss function (validation)
    """
    spec_name = spec_dict["name"]
    model_inputs = []
    model_obs = []

    # The first step is to generate an observable layer for each given datasets
    for dataset_dict in spec_dict["datasets"]:
        dataname = dataset_dict["name"]

        # Choose which kind of observable are we dealing with, since this defines the convolution
        if dataset_dict["hadronic"]:
            Obs_Layer = DY
        else:
            Obs_Layer = DIS

        # Define the operation (if any)
        op = operations.c_to_py_fun(dataset_dict["operation"], name=dataname)
        # Retrieve the list of fktables
        fktable_list = dataset_dict["fktables"]

        obs_list = []
        input_list = []
        # Now generate an input and output layer for each sub_obs of the dataset
        for i, fktable_dict in enumerate(fktable_list):
            # Input layer
            input_layer = operations.numpy_to_input(fktable_dict["xgrid"].T)
            input_list.append(input_layer)
            # Output layer
            obs_layer = Obs_Layer(
                output_dim=fktable_dict["ndata"],
                fktable=fktable_dict["fktable"],
                basis=fktable_dict["basis"],
                name="{0}_{1}".format(dataname, i),
                input_shape=(14,),
            )
            obs_list.append((input_layer, obs_layer))

        # Add the inputs to the lists of inputs of the model
        model_inputs += input_list
        # Append a combination of the operation to be applied (op) to the list
        # and the list of observable to which we want to applied the op
        model_obs.append((op, obs_list))

    def final_obs(pdf_layer):
        all_ops = []
        for operation, observables in model_obs:
            all_obs = []
            for i_layer, o_layer in observables:
                all_obs.append(o_layer(pdf_layer(i_layer)))
            all_ops.append(operation(all_obs))
        if len(all_ops) == 1:
            return all_ops[0]
        else:
            return concatenate(all_ops, axis=0, name=spec_name + "_full")

    if spec_dict["positivity"]:
        max_lambda = spec_dict["lambda"]
        if not positivity_initial:
            initial_lambda = max_lambda / pow(positivity_multiplier, positivity_steps)
        else:
            initial_lambda = positivity_initial
        out_tr_mask = Mask(
            bool_mask=spec_dict["trmask"], c=initial_lambda, name=spec_name
        )

        def out_tr_positivity(pdf_layer):
            return out_tr_mask(final_obs(pdf_layer))

        layer_info = {
            "inputs": model_inputs,
            "output_tr": out_tr_positivity,
            "loss_tr": losses.l_positivity(),
        }
        return layer_info

    # Now generate the mask layers to be applied to training and validation
    out_tr_mask = Mask(bool_mask=spec_dict["trmask"], name=spec_name)
    out_vl_mask = Mask(bool_mask=spec_dict["vlmask"], name=spec_name + "_val")

    def out_tr(pdf_layer):
        return out_tr_mask(final_obs(pdf_layer))

    def out_vl(pdf_layer):
        return out_vl_mask(final_obs(pdf_layer))

    # Generate the loss function as usual
    invcovmat = spec_dict["invcovmat_true"]
    loss = losses.l_invcovmat(invcovmat)

    invcovmat_tr = spec_dict["invcovmat"]
    loss_tr = losses.l_invcovmat(invcovmat_tr)

    invcovmat_vl = spec_dict["invcovmat_vl"]
    loss_vl = losses.l_invcovmat(invcovmat_vl)

    layer_info = {
        "inputs": model_inputs,
        "output": final_obs,
        "loss": loss,
        "output_tr": out_tr,
        "loss_tr": loss_tr,
        "output_vl": out_vl,
        "loss_vl": loss_vl,
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
            initializers.append(
                MetaLayer.select_initializer(initializer_name, seed=current_seed)
            )
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
            concatenator = base_layer_selector("concatenate")

            def output_layer(ilayer):
                result = layer(ilayer)
                return concatenator(result)

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
    out=14,
    seed=None,
    dropout=0.0,
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
        `inp`
            dimension of the xgrid. If inp=2, turns the x point into a (x, log(x)) pair
        `nodes'
            list of the number of nodes per layer of the PDF NN. Default: [15,8]
        `activation`
            list of activation functions to apply to each layer. Default: ["tanh", "linear"]
            if the number of activation function does not match the number of layers, it will add
            copies of the first activation function found
        `initializer_name`
            selects the initializer of the weights of the NN. Default: glorot_normal
        `layer_type`
            selects the type of architecture of the NN. Default: dense
        `flav_info`
            dictionary containing the information about each PDF (basis dictionary in the runcard)
            to be used by Preprocessing
        `out`
            number of output flavours of the model
        `seed`
            seed to initialize the NN
        `dropout`
            rate of dropout layer by layer

    Returns
    -------
        `layer_pdf`
            a function which, upon calling it with a tensor,
            will connect all PDF layers and output a tensor of size (batch_size, `out`)
        `dict_layers`
            a dictionary containing some of the intermediate layers
            (necessary for debugging and to compute intermediate quantities)
    """
    if nodes is None:
        nodes = [15, 8]
    if activations is None:
        activations = ["tanh", "linear"]
    # Safety check
    number_of_layers = len(nodes)
    number_of_activations = len(activations)
    if number_of_layers != number_of_activations:
        raise ValueError(
            "Number of activation functions does not match number of layers @ model_gen.py"
        )
    last_layer_nodes = nodes[-1]

    if layer_type == "dense":
        list_of_pdf_layers = generate_dense_network(
            inp, nodes, activations, initializer_name, seed=seed, dropout_rate=dropout
        )
    elif layer_type == "dense_per_flavour":
        # Define the basis size attending to the last layer in the network
        # TODO: this information should come from the basis information
        #       once the basis information is passed to this class
        list_of_pdf_layers = generate_dense_per_flavour_network(
            inp,
            nodes,
            activations,
            initializer_name,
            seed=seed,
            basis_size=last_layer_nodes,
        )

    # If the input is of type (x, logx)
    # create a x --> (x, logx) layer to preppend to everything
    if inp == 2:
        add_log = Lambda(lambda x: concatenate([x, operations.op_log(x)], axis=-1))

    def dense_me(x):
        """ Takes an input tensor `x` and applies all layers
        from the `list_of_pdf_layers` in order """
        if inp == 1:
            curr_fun = list_of_pdf_layers[0](x)
        else:
            curr_fun = list_of_pdf_layers[0](add_log(x))

        for dense_layer in list_of_pdf_layers[1:]:
            curr_fun = dense_layer(curr_fun)
        return curr_fun

    # Preprocessing layer (will be multiplied to the last of the denses)
    preproseed = seed + number_of_layers
    layer_preproc = Preprocessing(
        input_shape=(1,), name="pdf_prepro", flav_info=flav_info, seed=preproseed
    )

    # Apply preprocessing
    def layer_fitbasis(x):
        return operations.op_multiply([dense_me(x), layer_preproc(x)])

    # Evolution layer
    layer_evln = Rotation(input_shape=(last_layer_nodes,), output_dim=out)

    # Rotation layer, changes from the 8-basis to the 14-basis
    def layer_pdf(x):
        return layer_evln(layer_fitbasis(x))

    dict_layers = {
        "denses": dense_me,  # The set of the N dense layers
        "preprocessing": layer_preproc,  # The layer that applies preprocessing
        "fitbasis": layer_fitbasis,  # Applied preprocessing to the output of the denses
    }

    return layer_pdf, dict_layers
