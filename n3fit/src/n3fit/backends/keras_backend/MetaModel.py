"""
    MetaModel class

    Extension of the backend Model class containing some wrappers in order to absorb other
    backend-dependent calls.
"""

import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras import optimizers as Kopt
from n3fit.backends.keras_backend.operations import numpy_to_tensor

# Check the TF version to check if legacy-mode is needed (TF < 2.2)
tf_version = tf.__version__.split('.')
if int(tf_version[0]) == 2 and int(tf_version[1]) < 2:
    LEGACY = True
else:
    LEGACY = False

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
}

# Some keys need to work for everyone
for k, v in optimizers.items():
    v[1]["clipnorm"] = 1.0


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
                self.x_in[name] = numpy_to_tensor(input_tensor.tensor_content)
                self.tensors_in[name] = input_tensor
            except AttributeError:
                self.x_in[name] = None
                self.tensors_in[name] = None

        self.all_inputs = input_list
        self.all_outputs = output_list
        self.target_tensors = None
        self.eval_fun = None

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
        history = super().fit(x=x, y=y, epochs=epochs, **kwargs,)
        loss_dict = history.history
        return loss_dict

    def predict(self, x=None, **kwargs):
        """ Call super().predict with the right input arguments """
        x = self._parse_input(x)
        result = super().predict(x=x, **kwargs)
        return result

    def compute_losses(self):
        """
        This function is the fast-equivalent to the model ``evaluate(x,y)`` method.

        On first call it calls ``.evaluate(return_dict=True, verbose=0)`` to force
        the initialization of the test function.
        Subsequent calls of this method will (when applicable)
        directly call the internal evaluation function ``eval_fun``.
        This bypasses the pre- and post- evaluation steps, resulting in a ~10% speed up
        with respect to ``.evaluate(...)``

        Returns
        -------
            dict
                a dictionary with all partial losses of the model
        """
        if self.eval_fun is None:
            # We still need to perform some initialization
            if LEGACY:
                # For TF < 2.2 we need to generate the test_function ourselves
                self.make_test_function()
            else:
                return self.evaluate(return_dict=True, verbose=False)
        if LEGACY:
            # For tF < 2.2 we need to force the output to be a float
            ret = self.eval_fun()
            ret['loss'] = ret['loss'].numpy()
            return ret
        else:
            return self.eval_fun()

    def evaluate(self, x=None, y=None, **kwargs):
        """
        Wrapper around evaluate to take into account the case in which the data is already known
        when the model is compiled.
        """
        x = self._parse_input(self.x_in)
        if LEGACY and y is None:
            y = self.target_tensors
        result = super().evaluate(x=x, y=y, **kwargs)
        return result

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

        opt_function = opt_tuple[0]
        opt_args = opt_tuple[1]

        user_selected_args = {"learning_rate": learning_rate, "clipnorm": clipnorm}

        # Override defaults with user provided values
        for key, value in user_selected_args.items():
            if key in opt_args.keys() and value is not None:
                opt_args[key] = value

        # Instantiate the optimizer
        opt = opt_function(**opt_args)

        # If given target output, compile it together with the model for better performance
        if target_output is not None:
            if not isinstance(target_output, list):
                target_output = [target_output]
            # Tensorize
            self.target_tensors = target_output

        # Reset the evaluation function (if any)
        self.eval_fun = None

        super(MetaModel, self).compile(optimizer=opt, loss=loss)

    def make_test_function(self):
        """ If the model has been compiled with target data, it creates
        a specific evaluate function with the target data already evaluated.
        Otherwise return the normal tensorflow behaviour.
        """
        if self.eval_fun is not None:
            return self.eval_fun

        if self.target_tensors is None:
            return super().make_test_function()

        # Recover the target tensors and their lengths, we cannot rely
        # directly on the output from the model as we might have target_tensors
        # with 0 data points (if the tr/vl mask covers the whole set)
        lens = []
        tt = []
        for target in self.target_tensors:
            lens.append(target.size)
            tt.append(numpy_to_tensor(target))
        # Save target_tensors as tensors, as it might be useful for LEGACY
        self.target_tensors = tt

        # Get the name of the output layer
        # and add the suffix _loss to match TF behaviour
        out_names = [f"{i}_loss" for i in self.output_names]
        out_names.insert(0, "loss")

        @tf.function
        def eval_fun(*args):
            predictions = self(self._parse_input(None))
            # Concatenate the output to split them again as a list
            ypred = tf.concat(predictions, axis=-1)
            predspl = tf.split(ypred, lens, axis=-1)
            loss_list = [lfun(target, pred) for target, pred, lfun in zip(tt, predspl, self.loss)]
            ret = [tf.reduce_sum(loss_list)] + loss_list
            return dict(zip(out_names, ret))

        # Save the function so we don't go through this again
        self.eval_fun = eval_fun

        return eval_fun

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
