"""
    MetaModel class

    Extension of the backend Model class containing some wrappers in order to absorb other
    backend-dependent calls
"""
from keras.models import Model
import keras.optimizers as Kopt

from n3fit.backends.keras_backend.operations import numpy_to_tensor, numpy_to_input


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

    def fit(self, epochs = 1, **kwargs):
        if self.data_set:
            return super().fit(x = None, y = None, steps_per_epoch = 1, batch_size = None,  epochs = epochs,  **kwargs )
        else:
            return super().fit(epochs = epochs, **kwargs )
    def evaluate(self, **kwargs):
        if self.data_set:
            return super().evaluate(x = None, y = None, steps = 1, **kwargs)
        else:
            return super().evaluate()

    def __init__(self, input_tensors, output_tensors, extra_tensors=None, **kwargs):
        """
        This class behaves as keras.models.Model with the add on of extra_tensors.

        It is assumed here and elsewhere that the input_tensors are indeed tensors and
        not placeholders.
        The usage of placeholders will work as long as nothing outside of the ordinary is done
        i.e., it will work in that this class inherits from keras Model but all the add-ons
        of this class assume the inputs are tensors.
        
        The outputs however can be freely fixed (useful when fitting to known data)

        if extra_tensors are given, they should come in a list of tuples such that: 
            (input: np.array, output_layer: function)
        so that they can be fed to keras as output_layer(*input)
        The inputs of the extra_tensors will be turned into keras inputs.
        """
        self.data_set = None

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
                        inputs.append(numpy_to_input(i))
                else:
                    inputs = [numpy_to_input(ii)]
                o_tensor = oo(*inputs)
                input_list += inputs
                output_list.append(o_tensor)

        super(MetaModel, self).__init__(input_list, output_list, **kwargs)

    def compile(self, optimizer_name="RMSprop", learning_rate=0.05, loss=None, target_output = None, **kwargs):
        """
        Compile the model given an optimizer and a list of loss functions.
        The optimizer must be one of those implemented in the `optimizer` attribute of this class.

        Optionally a learning rate and a list of target outpout can be defined.

        # Arguments:
            - `optimizer_name`: string defining the optimizer to be used
            - `learning_rate`: learning rate of of the optimizer (if accepted as an argument, if not it will be ignored)
            - `loss` : list of loss functions to be pass to the model
            - `target_output`: list of outputs to compare the results to during fitting/evaluation
                               if given further calls to fit/evaluate must be done with y = None.
        """
        try:
            opt_tuple = self.optimizers[optimizer_name]
        except KeyError as e:
            raise NotImplementedError(
                f"[MetaModel.select_initializer] optimizer not implemented: {optimizer_name}"
            ) from e

        opt_function = opt_tuple[0]
        opt_args = opt_tuple[1]

        if "lr" in opt_args.keys():
            opt_args["lr"] = learning_rate

        opt_args["clipnorm"] = 1.0
        opt = opt_function(**opt_args)

        if target_output is not None:
            self.data_set = True

        super(MetaModel, self).compile(optimizer=opt, loss=loss, target_tensors = target_output,**kwargs)
