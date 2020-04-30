"""
    MetaModel class

    Extension of the backend Model class containing some wrappers in order to absorb other
    backend-dependent calls
"""

from tensorflow.keras.models import Model
from tensorflow.keras import optimizers as Kopt
from n3fit.backends.keras_backend import operations as op

# Define in this dictionary new optimizers as well as the arguments they accept
# (with default values if needed be)
optimizers = {
    "RMSprop": (Kopt.RMSprop, {"lr": 0.01}),
    "Adam": (Kopt.Adam, {"lr": 0.01}),
    "Adagrad": (Kopt.Adagrad, {}),
    "Adadelta": (Kopt.Adadelta, {"lr": 1.0}),
    "Adamax": (Kopt.Adamax, {}),
    "Nadam": (Kopt.Nadam, {}),
    "Amsgrad": (Kopt.Adam, {"lr": 0.01, "amsgrad": True}),
}


def _fill_placeholders(original_input, new_input=None):
    """
    Fills the placeholders of the original input with a new set of input

    Parameters
    ----------
        original_input: dictionary
            dictionary of input layers, can contain None
        new_input: list or dictionary
            list or dictionary of layers to substitute the None with
    """
    if new_input is None:
        return original_input
    x = {}
    i = 0
    for key, value in original_input.items():
        if value is None:
            try:
                x[key] = new_input[key]
            except TypeError:
                x[key] = new_input[i]
                i += 1
        else:
            x[key] = value
    return x


class MetaModel(Model):
    """
    The `MetaModel` behaves as the tensorflow.keras.model.Model class,
    with the addition of `tensor_content`:

    - tensor_content:
    Sometimes when fitting a network the input is fixed, in this case the input can be given
    together with the input_tensors by setting a `tensor_content` equal to the input value.
    This is done automatically when using the `numpy_to_input` function from
    `n3fit.backends.keras_backend.operations`

    Parameters
    ----------
        input_tensors: tensorflow.keras.layers.Input
            Input layer
        output_tensors: tensorflow.keras.layers.Layer
            Output layer
        **kwargs:
            keyword arguments to pass directly to Model
    """

    def __init__(self, input_tensors, output_tensors, **kwargs):
        self.has_dataset = False

        input_list = input_tensors
        output_list = output_tensors

        if not isinstance(input_list, list):
            input_list = [input_list]
        if not isinstance(output_list, list):
            output_list = [output_list]

        super(MetaModel, self).__init__(input_list, output_list, **kwargs)
        self.x_in = {}
        self.tensors_in = {}
        for input_tensor in input_list:
            # If the input contains a tensor_content, store it to use at predict/fit/eval times
            # otherwise, put a placeholder None as it will come from the outside
            name = input_tensor.op.name
            try:
                self.x_in[name] = input_tensor.tensor_content
                self.tensors_in[name] = input_tensor
            except AttributeError:
                self.x_in[name] = None
                self.tensors_in[name] = None

        self.all_inputs = input_list
        self.all_outputs = output_list
        self.target_tensors = None

    def _parse_input(self, extra_input, pass_numpy=True):
        """ Introduces the extra_input in the places asigned to the
        placeholders """
        if pass_numpy:
            return _fill_placeholders(self.x_in, extra_input)
        else:
            return _fill_placeholders(self.tensors_in, extra_input)

    def reinitialize(self):
        """ Run through all layers and reinitialize the ones that can be reinitialied """
        for layer in self.layers:
            if hasattr(layer, "reinitialize"):
                layer.reinitialize()

    def perform_fit(self, x=None, y=None, epochs=1, **kwargs):
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
        x = self._parse_input(x)
        if y is None:
            y = self.target_tensors
        history = super().fit(x=x, y=y, epochs=epochs, **kwargs,)
        loss_dict = history.history
        return loss_dict

    def predict(self, x=None, **kwargs):
        """ Call super().predict with the right input arguments """
        x = self._parse_input(x)
        result = super().predict(x=x, **kwargs)
        return result

    def compute_losses(self, *args, **kwargs):
        """
        Performs keras.evaluate and returns a dictionary containing the loss function for
        each of the outputs.

        This function acts as wrapper around Keras' evaluate and uses the information from
        the Keras' metric in order to return information on the losses in an unified format.

        Instead, the Keras' evaluate method returns either a list (with no information about which
        loss corresponds to which output) or a float (if there is a single output).

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
        x = self._parse_input(x)
        if y is None:
            y = self.target_tensors
            # TODO Ensure that no x or y were passed
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
            opt_tuple = optimizers[optimizer_name]
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
            if not isinstance(target_output, list):
                target_output = [target_output]
            self.target_tensors = None # TODO TF 2.2 target_output
            # Tensorize the input
            target = [op.numpy_to_tensor(i) for i in target_output]
        super(MetaModel, self).compile(optimizer=opt, target_tensors=target, loss=loss)

    def set_masks_to(self, names, val=0.0):
        """ Set all mask value to the selected value
        Masks in MetaModel should be named {name}_mask

        Mask are layers with one single weight (shape=(1,)) that multiplies the input

        Parameters
        ----------
            names: list
                list of masks to look for
            val: float
                selected value of the mask
        """
        mask_val = [val]
        for name in names:
            mask_name = f"{name}_mask"
            mask_w = self.get_layer(mask_name).weights[0]
            mask_w.assign(mask_val)

    def multiply_weights(self, layer_names, multiplier):
        """ Multiply all weights for the given layers by some scalar

        Parameters
        ----------
            layer_names: list
                list of names of the layers to update weights
            multiplier: float
                scalar number to multiply the weights by
        """
        for layer_name in layer_names:
            layer = self.get_layer(layer_name)
            w_val = layer.get_weights()
            w_ref = layer.weights
            for val, tensor in zip(w_val, w_ref):
                tensor.assign(val * multiplier)

    def apply_as_layer(self, x):
        """ Apply the model as a layer """
        x = self._parse_input(x, pass_numpy=False)
        return super().__call__(x)
