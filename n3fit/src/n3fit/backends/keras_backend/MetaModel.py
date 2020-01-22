"""
    MetaModel class

    Extension of the backend Model class containing some wrappers in order to absorb other
    backend-dependent calls
"""
from keras.models import Model, Sequential
import keras.optimizers as Kopt

from n3fit.backends.keras_backend.operations import numpy_to_input


class MetaModel(Model):
    """
    The goal of this class is to absorb all keras dependent code
    """

    # Define in this dictionary new optimizers as well as the arguments they accept
    # (with default values if needed be)
    optimizers = {
        "RMSprop": (Kopt.RMSprop, {"lr": 0.01}),
        "Adam": (Kopt.Adam, {"lr": 0.01}),
        "Adagrad": (Kopt.Adagrad, {}),
        "Adadelta": (Kopt.Adadelta, {}),
        "Adamax": (Kopt.Adamax, {}),
        "Nadam": (Kopt.Nadam, {}),
        "Amsgrad": (Kopt.Adam, {"lr": 0.01, "amsgrad": True}),
    }

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
        self.has_dataset = False

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

        self.all_inputs = input_list
        self.all_outputs = output_list
        self.internal_models = {}

        super(MetaModel, self).__init__(input_list, output_list, **kwargs)

    def perform_fit(self, epochs=1, **kwargs):
        """
        Performs forward (and backwards) propagation for the model for a given number of epochs.
        If the model was compiled with input and output data, they will not be passed through

        Returns
        -------
            `loss_dict`: dict
                a dictionary with all partial losses of the model
        """
        if self.has_dataset:
            history = super().fit(
                x=None, y=None, steps_per_epoch=1, batch_size=None, epochs=epochs, **kwargs
            )
        else:
            history = super().fit(epochs=epochs, **kwargs)
        loss_dict = history.history
        return loss_dict

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

    def evaluate(self, x = None, y = None, **kwargs):
        """
        Wrapper around evaluate to take into account the case in which the data is already known
        at the time of `.compile`.
        In this case the number of steps must be always specified and the input of x and y must
        be set to `None`.
        """
        if self.has_dataset:
            # Ensure that no x or y were passed
            result = super().evaluate(x=None, y=None, steps=1, **kwargs)
        else:
            result = super().evaluate(x = x, y = y, **kwargs)
        return result

    def compile(
        self, optimizer_name="RMSprop", learning_rate=0.05, loss=None, target_output=None, **kwargs
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

        if isinstance(opt_function, str):
	    # This allows for quickly drawing new optimizers that Keras might implement
            opt = opt_function
        else:
            opt = opt_function(**opt_args)

        if target_output is not None:
            self.has_dataset = True

        super(MetaModel, self).compile(
            optimizer=opt, loss=loss, target_tensors=target_output, **kwargs
        )

    def multiply_weights(self, key, layer_names, multiplier):
        """ Multiply all weights for the given layers by some scalar

        Parameters
        ----------
            `layer_names`: list
                list of names of the layers to update weights
            `multiplier`: float
                scalar number to multiply the weights by
        """
        internal_model = self.internal_models.get(key)
        if not internal_model:
            # Create an internal model to access the weights of the
            # layer we want to update
            # is this a bug or a feature?
            layers = [self.get_layer(i) for i in layer_names]
            internal_model = Sequential(layers)
            self.internal_models[key] = internal_model
        current_weights = internal_model.get_weights()
        new_weights = [i*multiplier for i in current_weights]
        internal_model.set_weights(new_weights)
