"""
    For a layer to be used by n3fit it should be contained in the layers dictionary below
    This dictionary has the following structure:

    'name of the layer' : ( Layer_class, {dictionary of arguments: defaults} )
"""

from keras.layers import Dense, Lambda, LSTM, Dropout, concatenate
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

def dense_per_flavour(basis_size = 8, concatenate_now = False, **dense_kwargs):
    """
    """

    # Need to generate a Dense layer per element of basis_size
    dense_basis = [base_layer_selector("dense", **dense_kwargs) for i in range(basis_size)]

    def apply_dense(xinput):
        """
        The input can be one single layer of a list of layer.
        If a single layer is given, all denses will "eat" it
        If, instead, a list of layers is given it should a) be of `basis_size`
        b) each dense of the basis will apply to only one of the layers
        """
        if isinstance(xinput, list): # TODO I'm not sure whether tensors will ever look as a duck
            if len(xinput) != basis_size:
                raise ValueError("You are evil")
            results = []
            for input_layer, den in zip(xinput, dense_basis):
                results.append(den(input_layer))
        else:
            results = [den(xinput) for den in dense_basis]

        # If we are in the last one, we need to concatenate
        if concatenate_now:
            output_layer = concatenate(results)
            return output_layer
        else:
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
    "dense_per_flavour" : (
        dense_per_flavour,
        {
            "input_shape": (1,),
            "kernel_initializer": "glorot_normal",
            "units": 5,
            "activation": "sigmoid",
            "basis_size" : 8,
            "concatenate_now" : False,
        },
    ),
    "LSTM": (
        LSTM_modified,
        {"kernel_initializer": "glorot_normal", "units": 5, "activation": "sigmoid"},
    ),
    "dropout": (Dropout, {"rate": 0.0}),
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
            "Layer not implemented in keras_backend/base_layers.py: {0}".format(layer_name)
        ) from e

    layer_class = layer_tuple[0]
    layer_args = layer_tuple[1]

    for key, value in kwargs.items():
        if key in layer_args.keys():
            layer_args[key] = value

    return layer_class(**layer_args)
