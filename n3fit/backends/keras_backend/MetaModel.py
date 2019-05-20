"""
    MetaModel class

    Extension of the backend Model class containing some wrappers in order to absorb other
    backend-dependent calls
"""
import numpy as np
from keras.models import Model
from keras.layers import Input
import keras.optimizers as Kopt
import keras.backend as K


class MetaModel(Model):
    """
    The goal of this class is to absorb all keras dependent code
    """

    # Define in this dictionary new optimizers as well as the arguments they accept (with default values if needed be)
    optimizers = {
            'RMSprop'  : (
                Kopt.RMSprop, {'lr': 0.01}
                ),
            'Adam'     : (
                Kopt.Adam, {'lr' : 0.01}
                ),
            'Adagrad'  : (
                Kopt.Adagrad, {} ),
            'Adadelta' : (
                Kopt.Adadelta, {} ),
            'Adamax'   : (
                Kopt.Adamax, {} ),
            'Nadam'    : (
                Kopt.Nadam, {} ),
            'Amsgrad'  : (
                Kopt.Adam, {'lr' : 0.01, 'amsgrad' : True }
                ),
           }

    def __init__(self, input_tensors, output_tensors, extra_tensors = None):
        """
        This class behaves as keras.models.Model with some add-ons:

        if extra_tensors are given, they should come in the combination
        (input: np.array, output_layer: function) so that they can be fed to keras as output_layer(*input)

        This function will turn all inputs into keras tensors
        """
        input_list = input_tensors
        output_list = output_tensors

        if not isinstance(input_list, list):
            input_list = [input_list]
        if not isinstance(output_list, list):
            output_list = [output_list]

        # Add extra tensors
        if extra_tensors is not None:
            for ii, oo in extra_tensors:
                inputs = []
                if isinstance(ii, list):
                    for i in ii:
                        inputs.append( self._np_to_input(i) )
                else:
                    inputs = [ self._np_to_input(ii) ]
                o_tensor = oo(*inputs)
                input_list += inputs
                output_list.append(o_tensor)

        super (MetaModel, self).__init__( input_list, output_list )

    def _np_to_input(self, x):
        """
        If x is a numpy array, make it into a numpy tensor
        """
        if isinstance(x, np.ndarray):
            tensor = K.constant(x)
            return Input( tensor = tensor )
        else:
            return x

    def compile(self, optimizer_name = 'RMSprop', learning_rate = 0.05, loss = None, **kwargs):
        """
        Compile the model given:
            - Optimizer
            - Learning Rate
            - List of loss functions
        """
        try:
            opt_tuple = self.optimizers[optimizer_name]
        except KeyError:
            raise Exception(f'[MetaModel.compile] Optimizer not found {optimizer_name}')

        opt_function = opt_tuple[0]
        opt_args = opt_tuple[1]

        if 'lr' in opt_args.keys():
            opt_args['lr'] = learning_rate

        opt_args['clipnorm'] = 1.0
        opt = opt_function(**opt_args)

        super (MetaModel, self).compile(optimizer = opt, loss = loss, **kwargs)
