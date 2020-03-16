"""
    MetaModel class

    Extension of the backend Model class containing some wrappers in order to absorb other
    backend-dependent calls
"""

import tensorflow as tf
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras import optimizers as Kopt
from tensorflow.keras import backend as K

from n3fit.backends.keras_backend.operations import numpy_to_input, batchit

import numpy as np


class MetaModel(Model):
    """
    The `MetaModel` behaves as the tensorflow.keras.model.Model class, with some additions:

    1. tensor_content
    Sometimes when fitting a network the input is fixed, in this case the input can be given
    together with the input_tensors by setting a `tensor_content` equal to the input value.
    This is done automatically when using the `numpy_to_input` function from
    `n3fit.backends.keras_backend.operations`

    2. extra_tensors
    Generally, when instantiating a model: `Model(x,y)`, y is an output layer which is connected
    to x. Adding `extra_tensors = ([x1, x2], y1)` will change the input of `Model(x,y)` to
    `Model([x, x1, x2], [y, y1])` where y is connected with x and y1 is connected with (x1,x2)
    automatically.
    

    Parameters
    ----------
        input_tensors: tensorflow.keras.layers.Input
            Input layer
        output_tensors: tensorflow.keras.layers.Layer
            Output layer
        extra_tensors: tuple
            Tuple of ([inputs], output)
    """

    # Define in this dictionary new optimizers as well as the arguments they accept
    # (with default values if needed be)
    optimizers = {
        "RMSprop": (Kopt.RMSprop, {"lr": 0.01}),
        "Adam": (Kopt.Adam, {"lr": 0.01}),
        "Adagrad": (Kopt.Adagrad, {}),
        "Adadelta": (Kopt.Adadelta, {"lr":1.0}),
        "Adamax": (Kopt.Adamax, {}),
        "Nadam": (Kopt.Nadam, {}),
        "Amsgrad": (Kopt.Adam, {"lr": 0.01, "amsgrad": True}),
    }

    def __init__(self, input_tensors, output_tensors, extra_tensors=None, **kwargs):
        self.has_dataset = False

        input_list = input_tensors
        output_list = output_tensors

        if not isinstance(input_list, list):
            input_list = [input_list]
        if not isinstance(output_list, list):
            output_list = [output_list]

        # Add extra tensors
        if extra_tensors is not None:
            # Check whether we are using the original shape
            # or forcing batch one
            if input_list and hasattr(input_list[0], 'original_shape'):
                keep_shape = input_list[0].original_shape
            else:
                keep_shape = True
            for ii, oo in extra_tensors:
                inputs = []
                if isinstance(ii, list):
                    for i in ii:
                        inputs.append(numpy_to_input(i, no_reshape = keep_shape))
                else:
                    inputs = [numpy_to_input(ii)]
                output_layer = oo(*inputs)
                # If we are not keeping the original shape (i.e., we added a batch dimension)
                # add it also to the output layer
                if not keep_shape:
                    output_layer = batchit(output_layer)
                input_list += inputs
                output_list.append(output_layer)

        super(MetaModel, self).__init__(input_list, output_list, **kwargs)
        if hasattr(input_list[0], 'tensor_content'):
            self.x_in = [i.tensor_content for i in input_list]
        else:
            self.x_in = None
        self.all_inputs = input_list
        self.all_outputs = output_list

    def reinitialize(self):
        """ Run through all layers and reinitialize the ones that can be reinitialied """
        for layer in self.layers:
            if hasattr(layer, "reinitialize"):
                layer.reinitialize()


    def perform_fit(self, x=None, y=None, steps_per_epoch=1, **kwargs):
        """
        Performs forward (and backwards) propagation for the model for a given number of epochs.

        The output of this function consists on a dictionary that maps the names of the metrics
        of the model (the loss functions) to the partial losses.

        If the model was compiled with input and output data, they will not be passed through.
        In this case by default the number of `steps_per_epoch` will be set to 1

        ex:
            {'loss': [100], 'dataset_a_loss1' : [67], 'dataset_2_loss': [33]}

        Returns
        -------
            `loss_dict`: dict
                a dictionary with all partial losses of the model
        """
        if x is None:
            x = self.x_in
        if self.has_dataset:
            history = super().fit(
                x=x, y=None, batch_size=1, **kwargs,
            )
        else:
            history = super().fit(x=x, y=y, steps_per_epoch=steps_per_epoch, **kwargs)
        loss_dict = history.history
        return loss_dict
    
    def predict(self, x = None, *args, **kwargs):
        if x is None:
            x = self.x_in
        result = super().predict(x = x, *args, **kwargs)
        return result

    def compute_losses(self, *args, **kwargs):
        """
        Performs keras.evaluate and returns a dictionary containing the loss function for
        each of the outputs.

        This function acts as wrapper around Keras' evaluate and uses the information from
        the Keras' metric in order to return information on the losses in an unified format.

        Instead, the Keras' evaluate method returns either a list (with no information about which loss corresponds to
        which output) or a float (if there is a single output).

        This function is compatible with the dict output by the fit history object.

        Any parameters passed to this function will be directly passed down to Keras.

        Returns
        -------
            `loss_dict`: dict
                a dictionary with all partial losses of the model
        """
        result = self.evaluate(*args, **kwargs)
        # get the name of all losses
        metricas = self.metrics_names
        if isinstance(result, float):
            # if there is only one we have to game it
            result = [result, result]
            # get the name of the last layer before the loss
            last_name = self.layers[-1].name
            metricas.append(f"{last_name}_loss")

        # now make it into a dictionary so it looks like the history object
        loss_dict = dict(zip(metricas, result))
        return loss_dict

    def evaluate(self, x=None, y=None, **kwargs):
        """
        Wrapper around evaluate to take into account the case in which the data is already known
        at the time of `.compile`.
        In this case the number of steps must be always specified and the input of x and y must
        be set to `None`.
        """
        if x is None:
            x = self.x_in
        if self.has_dataset:
            # Ensure that no x or y were passed
            result = super().evaluate(x=x, y=None, **kwargs)
        else:
            result = super().evaluate(x=x, y=y, **kwargs)
        return result

    def compile(
        self,
        optimizer_name="RMSprop",
        learning_rate=0.05,
        loss=None,
        target_output=None,
        **kwargs,
    ):
        """
        Compile the model given an optimizer and a list of loss functions.
        The optimizer must be one of those implemented in the `optimizer` attribute of this class.

        Options:
            - A learning rate and a list of target outpout can be defined.
                These will be passed down to the optimizer.
            - A `target_output` can be defined. If done in this way
                (for instance because we know the target data will be the same for the whole fit)
                the data will be compiled together with the model and won't be necessary to
                input it again in .fit()

        Parameters
        ----------
            `optimizer_name`: str
                string defining the optimizer to be used
            `learning_rate`: float
                learning rate of of the optimizer
                (if accepted as an argument, if not it will be ignored)
            `loss` : list
                list of loss functions to be pass to the model
            `target_output`: list
                list of outputs to compare the results to during fitting/evaluation
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
            self.has_dataset = True

        if not isinstance(target_output, list):
            target_output = [target_output]

        super(MetaModel, self).compile(
            optimizer=opt, loss=loss, target_tensors=target_output, **kwargs
        )

    def multiply_weights(self, layer_names, multiplier):
        """ Multiply all weights for the given layers by some scalar

        Parameters
        ----------
            `layer_names`: list
                list of names of the layers to update weights
            `multiplier`: float
                scalar number to multiply the weights by
        """
        for layer_name in layer_names:
            layer = self.get_layer(layer_name)
            w_val = layer.get_weights()
            w_ref = layer.weights
            for val, tensor in zip(w_val, w_ref):
                tensor.assign(val * multiplier)
