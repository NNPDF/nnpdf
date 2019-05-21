"""
    Library of functions to generate the NN objects

    observable_generator:
        Generates the output layers
    pdfNN_layer_generator:
        Generates the PDF NN layer to be fitted
"""
from layers import DIS
from layers import DY
from layers import Mask
from layers import Preprocessing, Rotation

from backends import operations
from backends import losses
from backends import MetaLayer
from backends import base_layer_selector, concatenate, Lambda


def observable_generator(
    spec_dict, positivity_initial=None, positivity_multiplier=1.05, verbose=False, positivity_steps=300
):
    """
    Receives a spec object (can be either a dataset or an experiment)
    Takes some optional arguments:
        positivity_initial: if given, set this number as the positivity multiplier for epoch 1
        positivity_multiplier: how much the positivity increases every 100 steps
        positivity_steps: if positivity_initial is not given, computes the initial by assuming we want,
                          after 100**positivity_steps epochs, to have the lambda of the runcard
    returns a tuple with:
    (input_layer, obs_layer, loss function, exp_data)
        - input_layer: list of input layers
        - obs_layer: one single output layer
        - loss_function: one single loss function
        - exp_data: one single set of true data to compare to
    """
    spec_name = spec_dict["name"]
    model_inputs = []
    model_obs = []

    # The first step is to generate an observable layer for each given datasets
    for dataset_dict in spec_dict["datasets"]:
        dataname = dataset_dict["name"]

        # Choose which kind of observable are we dealing with, since this will define the convolution
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
        # Append the combination of observable and the operation to be applied to the list of model_obs
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
        out_tr_mask = Mask(bool_mask=spec_dict["trmask"], c=initial_lambda, name=spec_name)

        def out_tr(pdf_layer):
            return out_tr_mask(final_obs(pdf_layer))

        layer_info = {"inputs": model_inputs, "output_tr": out_tr, "loss_tr": losses.l_positivity()}
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
):
    """
    Generates a NN that takes as input the xgrid value and outputs the basis of 14 PDFs
    fk tables use
    """
    if nodes is None:
        nodes = [15, 8]
    if activations is None:
        activations = ["tanh", "linear"]
    # Safety check
    ln = len(nodes)
    la = len(activations)
    if ln > la:
        activations = (ln - la) * ["tanh"] + activations
        print("Did not received enough activation functions, appending copies of the first activation foudn")
    elif la > ln:
        activations = activations[-ln:]
        print("Received more activations than layers, using only the last ones")

    # If dropout == 0, don't add dropout layer
    if dropout == 0.0:
        dropme = -99
    else:  # Add the dropout to the second to last layer
        dropme = ln - 2

    # Layer generation
    dl = []
    pre = inp
    for i, (units, activation) in enumerate(zip(nodes, activations)):
        if i == dropme and i != 0:
            dl.append(base_layer_selector("dropout", rate=dropout))

        init = MetaLayer.select_initializer(initializer_name, seed=seed + i)

        arguments = {"kernel_initializer": init, "units": int(units), "activation": activation, "input_shape": (pre,)}
        # Arguments that are not used by a given layer are just dropped
        tmpl = base_layer_selector(layer_type, **arguments)

        dl.append(tmpl)
        pre = int(units)

    if inp == 2:
        add_log = Lambda(lambda x: concatenate([x, operations.op_log(x)], axis=1))

    def dense_me(x):
        if inp == 1:
            curr_fun = dl[0](x)
        else:
            curr_fun = dl[0](add_log(x))

        for dense_layer in dl[1:]:
            curr_fun = dense_layer(curr_fun)
        return curr_fun

    # Preprocessing layer (will be multiplied to the last of the denses)
    preproseed = seed + ln
    layer_preproc = Preprocessing(input_shape=(1,), name="pdf_prepro", flav_info=flav_info, seed=preproseed)

    # Evolution layer
    layer_evln = Rotation(input_shape=(pre,), output_dim=out)

    # Generate multiplier layer
    # multiplier = Multiply(name = "MultiplyPreproc")

    # Apply preprocessing
    def layer_fitbasis(x):
        return operations.op_multiply([dense_me(x), layer_preproc(x)])

    # Rotation layer, changes from the 8-basis to the 14-basis
    def layer_pdf(x):
        return layer_evln(layer_fitbasis(x))

    dict_layers = {
        "denses": dense_me,  # The set of the N dense layers
        "preprocessing": layer_preproc,  # The layer that applies preprocessing
        "fitbasis": layer_fitbasis,  # Applied preprocessing to the output of the denses
    }

    return layer_pdf, dict_layers
