"""
    MetaModel class

    Extension of the backend Model class containing some wrappers in order to absorb other
    backend-dependent calls.
"""

import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras import optimizers as Kopt
from tensorflow.python.keras.utils import tf_utils
from n3fit.backends.keras_backend import operations as op

# Define in this dictionary new optimizers as well as the arguments they accept
# (with default values if needed be)
optimizers = {
    "RMSprop": (Kopt.RMSprop, {"learning_rate": 0.01}),
    "Adam": (Kopt.Adam, {"learning_rate": 0.01}),
    "Adagrad": (Kopt.Adagrad, {}),
    "Adadelta": (Kopt.Adadelta, {"learning_rate": 1.0}),
    "Adamax": (Kopt.Adamax, {}),
    "Nadam": (Kopt.Nadam, {}),
    "Amsgrad": (Kopt.Adam, {"learning_rate": 0.01, "amsgrad": True}),
}

# Some keys need to work for everyone
for k, v in optimizers.items():
    v[1]["clipnorm"] = 1.0


def _default_loss(y_true, y_pred):
    """ Default loss to be used when the model is compiled with loss = Null
    (for instance if the prediction of the model is already the loss """
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

    def __init__(self, input_tensors, output_tensors, **kwargs):
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

        super(MetaModel, self).__init__(input_list, output_list, **kwargs)
        self.x_in = {}
        self.tensors_in = {}
        for input_tensor in input_list:
            # If the input contains a tensor_content, store it to use at predict/fit/eval times
            # otherwise, put a placeholder None as it will come from the outside
            name = input_tensor.name.rsplit(":",1)[0]
            try:
                self.x_in[name] = op.numpy_to_tensor(input_tensor.tensor_content)
                self.tensors_in[name] = input_tensor
            except AttributeError:
                self.x_in[name] = None
                self.tensors_in[name] = None

        self.all_inputs = input_list
        self.all_outputs = output_list
        self.target_tensors = None
        self.compute_losses_function = None

    def _parse_input(self, extra_input=None, pass_content=True):
        """ Returns the input tensors the model was compiled with.
        Introduces the extra_input in the places asigned to the
        placeholders.

        If ``pass_content`` is set to ``False``, pass the tensor object.
        """
        if pass_content:
            return _fill_placeholders(self.x_in, extra_input)
        else:
            return _fill_placeholders(self.tensors_in, extra_input)

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
        x = self._parse_input(self.x_in)
        if y is None:
            y = self.target_tensors
        history = self.fit(x=x, y=y, epochs=epochs, **kwargs)
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
        # TODO might not work for TF < 2.2
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

        # The output of this function is to be used by python (and numpy) so we need to convert
        # the tensorflow variable to python primitives or numpy arrays.
        # Undocumented TF function that converts all the tensors from the ret dictionary to numpy arrays
        # if it dissapears, equivalent for us to {k: i.numpy() for k, i in ret.items()}
        return tf_utils.to_numpy_or_python_type(ret)

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
            self.target_tensors = [np.zeros((1,1)) for i in self.output_shape]
        else:
            if not isinstance(target_output, list):
                target_output = [target_output]
            self.target_tensors = target_output

        super(MetaModel, self).compile(optimizer=opt, loss=loss)

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

    def reset_layer_weights_to(self, layer_names, reference_vals):
        """ Set weights for the given layer to the given reference values

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
        x = self._parse_input(x, pass_content=False)
        try:
            return super().__call__(x)
        except ValueError:
            # TF < 2.1
            # TF 2.0 seems to fail with ValueError when passing a dictionary as an input
            y = x.values()
            return super().__call__(y)
