"""
    For a layer to be used by n3fit it should be contained in the layers dictionary below
    This dictionary has the following structure:

    'name of the layer' : ( Layer_class, {dictionary of arguments: defaults} )
"""

from keras.layers import Dense, Lambda, LSTM, Dropout, Concatenate, concatenate
from keras.backend import expand_dims


def LSTM_modified(**kwargs):
    """
    LSTM asks for a sample X timestep X features kind of thing so we need to reshape the input
    """
    the_lstm = LSTM(**kwargs)
    ExpandDim = Lambda(lambda x: expand_dims(x, axis=-1))

    def ReshapedLSTM(input_tensor):
        if len(input_tensor.shape) == 2:
            reshaped = ExpandDim(input_tensor)
            return the_lstm(reshaped)
        else:
            return the_lstm(input_tensor)

    return ReshapedLSTM


def dense_per_flavour(basis_size=8, kernel_initializer="glorot_normal", **dense_kwargs):
    """
    Generates a list of layers which can take as an input either one single layer
    or a list of the same size
    If taking one single layer, this one single layer will be the input of every layer in the list.
    If taking a list of layer of the same size, each layer on the list will take
    as input the layer on the input list in the same position.

    Note that, if the initializer is seeded, it should be a list where the seed is different
    for each element.

    i.e., if `basis_size` is 3 and is taking as input one layer A the output will be:
        [B1(A), B2(A), B3(A)]
    if taking, instead, a list [A1, A2, A3] the output will be:
        [B1(A1), B2(A2), B3(A3)]
    """
    if isinstance(kernel_initializer, str):
        kernel_initializer = basis_size * [kernel_initializer]

    # Need to generate a list of dense layers
    dense_basis = [
        base_layer_selector("dense", kernel_initializer=initializer, **dense_kwargs)
        for initializer in kernel_initializer
    ]

    def apply_dense(xinput):
        """
        The input can be either one single layer or a list of layers of
        length `basis_size`

        If taking one single layer, this one single layer will be the input of every
        layer in the list.
        If taking a list of layer of the same size, each layer on the list will take
        as input the layer on the input list in the same position.
        """
        if isinstance(xinput, (list, tuple)):
            if len(xinput) != basis_size:
                raise ValueError(
                    f"""The input of the dense_per_flavour and the basis_size
doesn't match, got a list of length {len(xinput)} for a basis_size of {basis_size}"""
                )
            results = [dens(ilayer) for dens, ilayer in zip(dense_basis, xinput)]
        else:
            results = [dens(xinput) for dens in dense_basis]

        return results

    return apply_dense


layers = {
    "dense": (
        Dense,
        {
            "input_shape": (1,),
            "kernel_initializer": "glorot_normal",
            "units": 5,
            "activation": "sigmoid",
        },
    ),
    "dense_per_flavour": (
        dense_per_flavour,
        {
            "input_shape": (1,),
            "kernel_initializer": "glorot_normal",
            "units": 5,
            "activation": "sigmoid",
            "basis_size": 8,
        },
    ),
    "LSTM": (
        LSTM_modified,
        {"kernel_initializer": "glorot_normal", "units": 5, "activation": "sigmoid"},
    ),
    "dropout": (Dropout, {"rate": 0.0}),
    "concatenate": (Concatenate, {}),
}


def base_layer_selector(layer_name, **kwargs):
    """
        Given a layer name, looks for it in the `layers` dictionary and returns an instance.

        The layer dictionary defines a number of defaults
        but they can be overwritten/enhanced through kwargs

        Parameters
        ----------
            `layer_name
                str with the name of the layer
            `**kwargs`
                extra optional arguments to pass to the layer (beyond their defaults)
    """
    try:
        layer_tuple = layers[layer_name]
    except KeyError as e:
        raise NotImplementedError(
            "Layer not implemented in keras_backend/base_layers.py: {0}".format(
                layer_name
            )
        ) from e

    layer_class = layer_tuple[0]
    layer_args = layer_tuple[1]

    for key, value in kwargs.items():
        if key in layer_args.keys():
            layer_args[key] = value

    return layer_class(**layer_args)
