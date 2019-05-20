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
    ExpandDim = Lambda(
            lambda x: expand_dims(x, axis = -1)
        )
    def ReshapedLSTM(input_tensor):
        if len(input_tensor.shape) == 2:
            reshaped = ExpandDim(input_tensor)
            return the_lstm(reshaped)
        else:
            return the_lstm(input_tensor)
    return ReshapedLSTM


layers = {
        'dense' : (Dense, {
            'input_shape' : (1,),
            'kernel_initializer' : 'glorot_normal',
            'units' : 5,
            'activation' : 'sigmoid',
            }),
        'LSTM' : (LSTM_modified, {
            'kernel_initializer' : 'glorot_normal',
            'units' : 5,
            'activation' : 'sigmoid',
            }),
        'dropout' : (Dropout, {
            'rate' : 0.0,
            }),
        }



def base_layer_selector(layer_name, **kwargs):
    try:
        layer_tuple = layers[layer_name]
    except KeyError:
        raise Exception("Layer not implemented in keras_backend/base_layers.py: {0}".format(layer_name))

    layer_class = layer_tuple[0]
    layer_args = layer_tuple[1]

    for key, value in kwargs.items():
        if key in layer_args.keys():
            layer_args[key] = value

    return layer_class( **layer_args )
