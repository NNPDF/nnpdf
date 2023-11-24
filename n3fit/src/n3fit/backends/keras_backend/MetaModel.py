"""
    MetaModel class

    Extension of the backend Model class containing some wrappers in order to absorb other
    backend-dependent calls.
"""

import re

import h5py
import numpy as np
import tensorflow as tf
from tensorflow.keras import optimizers as Kopt
from tensorflow.keras.models import Model
from tensorflow.python.keras.utils import tf_utils  # pylint: disable=no-name-in-module

import n3fit.backends.keras_backend.operations as op

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


def _default_loss(y_true, y_pred):  # pylint: disable=unused-argument
    """Default loss to be used when the model is compiled with loss = Null
    (for instance if the prediction of the model is already the loss"""
    return op.sum(y_pred)


class MetaModel(Model):
    """
    The model wraps keras.Model and adds some custom behaviour. Most notably it
    allows supplying constant values for input arguments, which are used when
    training and making predictions with the model (note that constants need to
    be explicitly registered as inputs, see
    https://github.com/keras-team/keras/issues/11912). These inputs can be
    passed in the ``input_values`` parameter, or gathered from the
    ``tensor_content`` attribute of the ``input_tensors``, which is set
    automatically when using the ``numpy_to_input`` function from
    :py:mod:`n3fit.backends.keras_backend.operations`.

    Parameters
    ----------
        input_tensors: dict[Any, tensorflow.keras.layers.Input]
            Input layer
        output_tensors: tensorflow.keras.layers.Layer
            Output layer
        input_values: dict[Any, array_like]
            Constant values for the input layer, to be supplied when making
            predictions with the model.

        **kwargs:
            keyword arguments to pass directly to Model
    """

    accepted_optimizers = optimizers

    def __init__(self, input_tensors, output_tensors, scaler=None, input_values=None, **kwargs):
        self.has_dataset = False
        self.required_slots = set()

        if input_values is None:
            input_values = {}

        if not isinstance(input_tensors, dict):
            raise TypeError("Expecting input_tensors to be a dict")

        if not isinstance(input_values, dict):
            raise TypeError("Expecting input_values to be a dict or None")

        x_in = {}
        # Go over the inputs. If we can deduce a constant value, either because
        # it is set in input_values or because it has a tensor_content, we
        # store it. Otherwise we mark the input as required when making
        # predictions.
        for k, v in input_tensors.items():
            if k in input_values:
                x_in[k] = input_values[k]
            elif hasattr(v, "tensor_content"):
                x_in[k] = op.numpy_to_tensor(v.tensor_content)
            else:
                self.required_slots.add(k)
        super().__init__(input_tensors, output_tensors, **kwargs)

        self.x_in = x_in
        self.input_tensors = input_tensors

        self.target_tensors = None
        self.compute_losses_function = None
        self._scaler = scaler

    @tf.autograph.experimental.do_not_convert
    def _parse_input(self, extra_input=None):
        """Returns the input data the model was compiled with.
        Introduces the extra_input in the places asigned to the placeholders.

        If the model was generated with a scaler, the input will be scaled accordingly
        """
        if extra_input is None:
            if self.required_slots:
                raise ValueError(f"The following inputs must be provided: {self.required_slots}")
            return self.x_in

        if not isinstance(extra_input, dict):
            raise TypeError("extra_input should be a dict or None")

        if diff := (self.required_slots - extra_input.keys()):
            raise ValueError(f"The following inputs must be provided {diff}")

        if diff := (extra_input.keys() - (self.x_in.keys() | self.required_slots)):
            raise ValueError(f"The following inputs are unknown {diff}")

        if self._scaler is not None:
            extra_input = {k: self._scaler(i) for k, i in extra_input.items()}

        return {**self.x_in, **extra_input}

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
        x_params = self._parse_input(x)
        if y is None:
            y = self.target_tensors
        history = super().fit(x=x_params, y=y, epochs=epochs, **kwargs)
        loss_dict = history.history
        return loss_dict

    def predict(self, x=None, **kwargs):
        """Call super().predict with the right input arguments"""
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

        super().compile(optimizer=opt, loss=loss)

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
        """Apply the model as a layer"""
        all_input = {**self.input_tensors, **x}
        return all_input, super().__call__(all_input)

    def get_layer_re(self, regex):
        """Get all layers matching the given regular expression"""
        check = lambda x: re.match(regex, x.name)
        return list(filter(check, self.layers))

    def get_replica_weights(self, i_replica):
        """
        Get the weights of replica i_replica.

        This assumes that the only weights are in layers called
        ``NN_{i_replica}`` and ``preprocessing_factor_{i_replica}``


        Parameters
        ----------
            i_replica: int

        Returns
        -------
            dict
                dictionary with the weights of the replica
        """
        NN_weights = [
            tf.Variable(w, name=w.name) for w in self.get_layer(f"NN_{i_replica}").weights
        ]
        prepro_weights = [
            tf.Variable(w, name=w.name)
            for w in self.get_layer(f"preprocessing_factor_{i_replica}").weights
        ]
        weights = {"NN": NN_weights, "preprocessing_factor": prepro_weights}

        return weights

    def set_replica_weights(self, weights, i_replica=0):
        """
        Set the weights of replica i_replica.

        This assumes that the only weights are in layers called
        ``NN_{i_replica}`` and ``preprocessing_factor_{i_replica}``

        Parameters
        ----------
            weights: dict
                dictionary with the weights of the replica
            i_replica: int
                the replica number to set, defaulting to 0
        """
        self.get_layer(f"NN_{i_replica}").set_weights(weights["NN"])
        self.get_layer(f"preprocessing_factor_{i_replica}").set_weights(
            weights["preprocessing_factor"]
        )

    def split_replicas(self):
        """
        Split the single multi-replica model into a list of separate single replica models,
        maintaining the current state of the weights.

        This is a dirty implementation that takes all replicas and extracts the one of interest.

        Returns
        -------
            list
                list of single replica models
        """
        replicas = []
        for i_replica, replica in enumerate(self.single_replicas.replicas):
            replica.set_replica_weights(self.get_replica_weights(i_replica))

            # pick single photon
            if "add_photons" in self.layers:
                replica.get_layer("add_photons").set_photon(
                    self.get_layer("add_photons").get_photon(i_replica)
                )
            replicas.append(replica)

        return replicas

    def load_identical_replicas(self, model_file):
        """
        From a single replica model, load the same weights into all replicas.
        """
        weights = self._format_weights_from_file(model_file)

        num_replicas = self.output.shape[-1]
        for i_replica in range(num_replicas):
            self.set_replica_weights(weights, i_replica)

    def _format_weights_from_file(self, model_file):
        """Read weights from a .h5 file and format into a dictionary of tf.Variables"""
        weights = {}

        with h5py.File(model_file, 'r') as f:
            # look at layers of the form NN_i and take the lowest i
            i_replica = 0
            while f"NN_{i_replica}" not in f:
                i_replica += 1

            weights["NN"] = self._extract_weights(f[f"NN_{i_replica}"], "NN", i_replica)
            weights["preprocessing_factor"] = self._extract_weights(
                f[f"preprocessing_factor_{i_replica}"], "preprocessing_factor", i_replica
            )

        return weights

    def _extract_weights(self, h5_group, weights_key, i_replica):
        """Extract weights from a h5py group, turning them into Tensorflow variables"""
        weights = []

        def append_weights(name, node):
            if isinstance(node, h5py.Dataset):
                weight_name = node.name.split("/", 2)[-1]
                weight_name = weight_name.replace(f"NN_{i_replica}", f"NN_0")
                weight_name = weight_name.replace(
                    f"preprocessing_factor_{i_replica}", f"preprocessing_factor_0"
                )
                weights.append(tf.Variable(node[()], name=weight_name))

        h5_group.visititems(append_weights)

        # have to put them in the same order
        weights_ordered = []
        weights_model_order = [w.name for w in self.get_replica_weights(0)[weights_key]]
        for w in weights_model_order:
            for w_h5 in weights:
                if w_h5.name == w:
                    weights_ordered.append(w_h5)

        return weights_ordered

    def attach_single_replicas(self, single_replicas):
        self.single_replicas = SingleReplicaList(single_replicas)


class SingleReplicaList:
    """
    Wrapper around list of single replicas to prevent Keras from picking them up as layers
    of the multi-replica model
    """

    def __init__(self, replicas):
        self.replicas = replicas
