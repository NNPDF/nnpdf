"""
    MetaModel class

    Extension of the backend Model class containing some wrappers in order to absorb other
    backend-dependent calls.
"""

import re
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras import optimizers as Kopt
from tensorflow.python.keras.utils import tf_utils  # pylint: disable=no-name-in-module
import n3fit.backends.keras_backend.operations as op

alphas_opt = tf.keras.optimizers.Adam(learning_rate=0)

# Check the TF version to check if legacy-mode is needed (TF < 2.2)
tf_version = tf.__version__.split(".")
if int(tf_version[0]) == 2 and int(tf_version[1]) < 2:
    raise NotImplementedError("n3fit needs TF > 2.2 in order to work")


# We need a function to transform tensors to numpy/python primitives
# which is not part of the official TF interface and can change with the version
if hasattr(tf_utils, "to_numpy_or_python_type"):
    _to_numpy_or_python_type = tf_utils.to_numpy_or_python_type
elif hasattr(tf_utils, "sync_to_numpy_or_python_type"):  # from TF 2.5
    _to_numpy_or_python_type = tf_utils.sync_to_numpy_or_python_type
else:  # in case of disaster
    _to_numpy_or_python_type = lambda ret: {k: i.numpy() for k, i in ret.items()}


# Define in this dictionary new optimizers as well as the arguments they accept
# (with default values if needed be)
optimizers = {
    "RMSprop": (Kopt.RMSprop, {"learning_rate": 0.01}),
    "Adam": (Kopt.Adam, {"learning_rate": 0.01}),
    "Adagrad": (Kopt.Adagrad, {}),
    "Adadelta": (Kopt.Adadelta, {"learning_rate": 1.0}),
    "Adamax": (Kopt.Adamax, {}),
    "Nadam": (Kopt.Nadam, {"learning_rate": 0.001}),
    "Amsgrad": (Kopt.Adam, {"learning_rate": 0.01, "amsgrad": True}),
    "SGD": (Kopt.SGD, {"learning_rate": 0.01, "momentum": 0.0, "nesterov": False}),
}

# Some keys need to work for everyone
for k, v in optimizers.items():
    v[1]["clipnorm"] = 1.0


def _default_loss(y_true, y_pred): # pylint: disable=unused-argument
    """Default loss to be used when the model is compiled with loss = Null
    (for instance if the prediction of the model is already the loss"""
    return op.sum(y_pred)


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

    accepted_optimizers = optimizers

    def __init__(self, input_tensors, output_tensors, scaler=None, **kwargs):
        self.has_dataset = False

        input_list = input_tensors
        output_list = output_tensors

        if isinstance(input_list, dict):
            # if this is a dictionary, convert it to a list for now
            input_list = input_tensors.values()
        elif not isinstance(input_list, list):
            # if it is not a dict but also not a list, make it into a 1-element list and pray
            input_list = [input_list]

        if isinstance(output_list, dict):
            # if this is a dictionary, convert it to a list for now
            output_list = output_tensors.values()
        elif not isinstance(output_list, list):
            # if it is not a dict but also not a list, make it into a 1-element list and pray
            output_list = [output_list]

        # Note: there used to be two possible options when creating a model:
        # - Give placeholder tensors (for which the content will be given at run time)
        # - Give tensors with content* (for which the content is stored with the model
        # *this option was dropped at some point by TF, the code below keeps this behaviour
        # We will store within the model the following quantities:
        #   -> x_in: arrays containing the x-input to the model
        #   -> tensors_in: when the x-input is known at compile time, we store a reference to the tensor
        # We pass TensorFlow a dictionary {k: tensor} containing placeholders which will be automatically filled
        # whenever x_in/tensor_in is known at compile time

        x_in = {}
        tensors_in = {}
        input_dict = {}
        for input_tensor in input_list:
            # If the input contains a tensor_content, store it to use at predict/fit/eval times
            # otherwise, put a placeholder None as it will come from the outside
            name = input_tensor.name.rsplit(":", 1)[0]
            input_dict[name] = input_tensor
            try:
                x_in[name] = op.numpy_to_tensor(input_tensor.tensor_content)
                tensors_in[name] = input_tensor
            except AttributeError:
                x_in[name] = None
                tensors_in[name] = None

        super().__init__(input_dict, output_list, **kwargs)

        self.x_in = x_in
        self.tensors_in = tensors_in

        self.target_tensors = None
        self.compute_losses_function = None
        self._scaler = scaler

    def _parse_input(self, extra_input=None):
        """Returns the input data the model was compiled with.
        Introduces the extra_input in the places asigned to the placeholders.

        If the model was generated with a scaler, the input will be scaled accordingly
        """
        if isinstance(extra_input, dict) and self._scaler is not None:
            extra_input = {k: self._scaler(i) for k, i in extra_input.items()}
        elif isinstance(extra_input, (tuple, list)) and self._scaler is not None:
            extra_input = [self._scaler(i) for i in extra_input]
        return _fill_placeholders(self.x_in, extra_input)

    def perform_fit(self, x=None, y=None, epochs=1, **kwargs):
        """
        Performs forward (and backwards) propagation for the model for a given number of epochs.

        The output of this function consists on a dictionary that maps the names of the metrics
        of the model (the loss functions) to the partial losses.

        If the model was compiled with input and output data, they will not be passed through.
        In this case by default the number of `epochs` will be set to 1

        ex:
            {'loss': [100], 'dataset_a_loss1' : [67], 'dataset_2_loss': [33]}

        Returns
        -------
            loss_dict: dict
                a dictionary with all partial losses of the model
        """
        x = self._parse_input(x)
        if y is None:
            y = self.target_tensors
        history = super().fit(x=x, y=y, epochs=epochs, **kwargs)
        loss_dict = history.history
        return loss_dict

    def predict(self, x=None, **kwargs):
        """ Call super().predict with the right input arguments """
        x = self._parse_input(x)
        result = super().predict(x=x, **kwargs)
        return result

    def compute_losses(self):
        """
        This function is equivalent to the model ``evaluate(x,y)`` method of most TensorFlow models
        which return a dictionary of losses per output layer.
        The losses reported in the ``evaluate`` method for n3fit are, however, summed over replicas.
        Instead the loss we are interested in is usually the output of the model (i.e., predict)
        This function then generates a dict of partial losses of the model separated per replica.
        i.e., the output for experiment {'LHC_exp'} will be an array of Nrep elements.

        Returns
        -------
            dict
                a dictionary with all partial losses of the model
        """
        if self.compute_losses_function is None:
            # If it is the first time we are passing through, compile the function and save it
            out_names = [f"{i}_loss" for i in self.output_names]
            out_names.insert(0, "loss")

            # Compile a evaluation function
            @tf.function
            def losses_fun():
                predictions = self(self._parse_input(None))
                # If we only have one dataset the output changes
                if len(out_names) == 2:
                    predictions = [predictions]
                total_loss = tf.reduce_sum(predictions, axis=0)
                ret = [total_loss] + predictions
                return dict(zip(out_names, ret))

            self.compute_losses_function = losses_fun

        ret = self.compute_losses_function()

        # The output of this function is to be used by python (and numpy)
        # so we need to convert the tensors
        return _to_numpy_or_python_type(ret)

    def compile(
        self,
        optimizer_name="RMSprop",
        learning_rate=None,
        loss=None,
        target_output=None,
        clipnorm=None,
        **kwargs,
    ):
        """
        Compile the model given an optimizer and a list of loss functions.
        The optimizer must be one of those implemented in the `optimizer` attribute of this class.

        Options:
            - A learning rate and a list of target outpout can be defined.
                These will be passed down to the optimizer.
            - A ``target_output`` can be defined. If done in this way
                (for instance because we know the target data will be the same for the whole fit)
                the data will be compiled together with the model and won't be necessary to
                input it again when calling the ``perform_fit`` or ``compute_losses`` methods.

        Parameters
        ----------
            optimizer_name: str
               string defining the optimizer to be used
            learning_rate: float
               learning rate of of the optimizer
               (if accepted as an argument, if not it will be ignored)
            loss: list
               list of loss functions to be pass to the model
            target_output: list
                list of outputs to compare the results to during fitting/evaluation
                if given further calls to fit/evaluate must be done with y = None.
        """
        try:
            opt_tuple = optimizers[optimizer_name]
        except KeyError as e:
            raise NotImplementedError(
                f"[MetaModel.select_initializer] optimizer not implemented: {optimizer_name}"
            ) from e

        if loss is None:
            loss = _default_loss

        opt_function = opt_tuple[0]
        opt_args = opt_tuple[1]

        user_selected_args = {"learning_rate": learning_rate, "clipnorm": clipnorm}

        # Override defaults with user provided values
        for key, value in user_selected_args.items():
            if key in opt_args.keys() and value is not None:
                opt_args[key] = value

        # Instantiate the optimizer
        opt = opt_function(**opt_args)

        # If given target output is None, target_output is unnecesary, save just a zero per output
        if target_output is None:
            self.target_tensors = [np.zeros((1, 1)) for i in self.output_shape]
        else:
            if not isinstance(target_output, list):
                target_output = [target_output]
            self.target_tensors = target_output

        from tensorflow_addons.optimizers import MultiOptimizer
        non_alphas_layers = [i for i in self.layers if i.name is not "alphas_layer"]
        alphas_layers = [i for i in self.layers if i.name is "alphas_layer"]
        opts_and_layers = [(opt, non_alphas_layers), (alphas_opt, alphas_layers)]
        multi_optimizers = MultiOptimizer(opts_and_layers)

        super().compile(optimizer=multi_optimizers, loss=loss)

    def set_masks_to(self, names, val=0.0):
        """Set all mask value to the selected value
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

    def reset_layer_weights_to(self, layer_names, reference_vals):
        """Set weights for the given layer to the given reference values

        The ``reference_vals`` values list must be a list of the same size
        of ``layer_names`` and it must consist of numpy arrays that perfectly
        align to the reference layer weights.
        In the special case of 1-weight layers it admits a scalar as input.

        Parameters
        ----------
            layer_names: list
                list of names of the layers to update weights
            reference_vals: list(float) or list(arrays)
                list of scalar or arrays to assign to each layer
        """
        for layer_name, values in zip(layer_names, reference_vals):
            if np.isscalar(values):
                values = np.array([[values]])
            layer = self.get_layer(layer_name)
            all_w = layer.weights
            for w, v in zip(all_w, values):
                w.assign(v)

    def apply_as_layer(self, x):
        """ Apply the model as a layer """
        all_input = _fill_placeholders(self.tensors_in, x)
        return all_input, super().__call__(all_input)

    def get_layer_re(self, regex):
        """ Get all layers matching the given regular expression """
        check = lambda x: re.match(regex, x.name)
        return list(filter(check, self.layers))
