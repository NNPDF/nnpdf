from keras.layers import add as keras_add
from keras.layers import subtract as keras_subtract
from keras.layers import Lambda as keras_Lambda
from keras.layers import multiply as keras_multiply

from keras.layers import Input
from keras.models import Model
from keras import backend as K

def numpy_to_input(numpy_array):
    return Input( tensor = K.constant(numpy_array))

# All operations must accept as an input a list or tuple of keras layers
# with the same output shape
# and return a keras operation with, again, the same output shape

def c_to_py_fun(op_name, name, default = 'ADD'):
    """
    Map between the NNPDF operations and the operations defined in this file
    Any new backend must implement such a mapping
    """ # TODO: shouldn't this be outside of the backend folder then
        # surely...
    d = {
        'NULL' : op_null,
        'ADD' : op_add,
        'RATIO' : op_ratio,
        'ASY' : op_asy,
        'SMN' : op_smn,
            }
    if op_name not in d.keys():
        print("Operation name not recognised, defaulting to {0}".format(default))
        return d[default]

    def operation_fun(o_list):
        return d[op_name]( o_list, name = name )

    return operation_fun

def op_multiply(o_list, **kwargs):
    """
    Receives a list of layers of the same output size and multiply them element-wise
    """
    return keras_multiply(o_list, **kwargs)

def op_multiply_dim(o_list, **kwargs):
    """
    Bypass in order to multiply two layers with different output dimension
    for instance: (10000 x 14) * (14)
    as the normal keras multiply don't accept it (but somewhow it does accept it doing it like this)
    """
    if len(o_list) != 2:
        raise Exception("The number of observables is incorrect, operations.py:op_multiply_dim")
    create_operation = keras_Lambda(lambda inputs: inputs[0] * inputs[1])
    return create_operation(o_list)

def op_null(o_list, **kwargs):
    """
    Not a compound object, do nothing
    """
    return o_list[0]

def op_add(o_list, **kwargs):
    """
    Sum a list of layers with the same output dim
    """
    return keras_add(o_list, **kwargs)

def op_subtract(o_list, **kwargs):
    """
    Subtract all observables
    """
    return keras_subtract(o_list)

def op_log(o_tensor, **kwargs):
    """
    Computes the logarithm of the input
    """
    return K.log(o_tensor)


def op_ratio(o_list, **kwargs):
    """
    Take the ratio of two observables
    """
    if len(o_list) != 2:
        raise Exception("The number of observables is incorrect, operations.py:op_ratio")

    division_layer = keras_Lambda(lambda inputs: inputs[0] / inputs[1], **kwargs)
    return division_layer(o_list)

def op_asy(o_list, **kwargs):
    """
    Perform the asymmetry operation on two observables
    """
    if len(o_list) != 2:
        raise Exception("The number of observables is incorrect, operations.py:op_asy")

    subtraction = keras_subtract(o_list)
    addition = op_add(o_list)
    return op_ratio( [subtraction, addition], **kwargs )

def op_smn(o_list, **kwargs):
    """
    Normalised sum
    """
    if len(o_list) != 4:
        raise Exception("The number of observables is incorrect, operations.py:op_smn")
    numer = op_add( o_list[:2] )
    denom = op_add( o_list[2:] )
    return op_ratio( [numer, denom], **kwargs )
    
