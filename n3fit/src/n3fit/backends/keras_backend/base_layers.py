"""
This module defines custom base layers to be used by the n3fit
Neural Network.
These layers can use the keras standard set of activation function
or implement their own.

For a layer to be used by n3fit it should be contained in the `layers` dictionary defined below.
This dictionary has the following structure:

    'name of the layer' : ( Layer_class, {dictionary of arguments: defaults} )

In order to add custom activation functions, they must be added to
the `custom_activations` dictionary with the following structure:

    'name of the activation' : function

The names of the layer and the activation function are the ones to be used in the n3fit runcard.
"""

import numpy as np
import keras.backend as K
from keras import ops as Kops
from keras import random as krandom
import tensorflow as tf
import math

from keras.layers import Dense as KerasDense
from keras.layers import Dropout, Lambda, Layer
from keras.layers import Input  # pylint: disable=unused-import
from keras.layers import LSTM, Concatenate
from keras.regularizers import l1_l2

from . import operations as ops
from .MetaLayer import MetaLayer
from contextlib import contextmanager


# Custom activation functions
def square_activation(x):
    """Squares the input"""
    return x * x


def square_singlet(x):
    """Square the singlet sector
    Defined as the two first values of the NN"""
    singlet_squared = x[..., :2] ** 2
    return ops.concatenate([singlet_squared, x[..., 2:]], axis=-1)


def modified_tanh(x):
    """A non-saturating version of the tanh function"""
    return ops.absolute(x) * ops.tanh(x)


def leaky_relu(x):
    """Computes the Leaky ReLU activation function"""
    return ops.leaky_relu(x, alpha=0.2)


custom_activations = {
    "square": square_activation,
    "square_singlet": square_singlet,
    "leaky_relu": leaky_relu,
    "modified_tanh": modified_tanh,
}


def LSTM_modified(**kwargs):
    """
    LSTM asks for a sample X timestep X features kind of thing so we need to reshape the input
    """
    the_lstm = LSTM(**kwargs)
    ExpandDim = Lambda(lambda x: ops.expand_dims(x, axis=-1))

    def ReshapedLSTM(input_tensor):
        if len(input_tensor.shape) == 2:
            reshaped = ExpandDim(input_tensor)
            return the_lstm(reshaped)
        else:
            return the_lstm(input_tensor)

    return ReshapedLSTM


class VBDense(Layer):
    """
    Mean-field variational Bayesian dense layer for n3fit (backend-agnostic).

    Mirrors the clean `VBLinear` reference (three explicit forward paths + one
    analytic-KL helper), written with `keras.ops` / `keras.random` so it runs on
    any Keras-3 backend.
    """

    def __init__(
        self,
        out_features: int,
        in_features: int,
        prior_prec: float = 1.0,
        std_init: float = -9.0,
        map: bool = False,
        bayesian_bias: bool = False,
    ):
        super().__init__()
        self.output_dim = out_features
        self.input_dim = in_features
        self.map = map
        self.prior_prec = tf.cast(prior_prec, K.floatx())
        self.std_init = tf.cast(std_init, K.floatx())
        self.bayesian_bias = bayesian_bias
        self.lbound = -30 if K.floatx() == 'float64' else -20
        self.ubound = 11
        self.eps = 1e-12 if K.floatx() == 'float64' else 1e-8
        self.training = True
        # Frozen eval draws (plain attributes, NOT tracked weights).
        self.random = None
        self.random_b = None

    def build(self, input_shape):
        self.bias = self.add_weight(
            name='bias',
            shape=(self.output_dim,),
            initializer='glorot_normal',
            trainable=True,
            dtype=K.floatx(),
        )
        self.mu_w = self.add_weight(
            name='mu_w',
            shape=(self.output_dim, self.input_dim),
            initializer='glorot_normal',
            trainable=True,
            dtype=K.floatx(),
        )
        self.logsig2_w = self.add_weight(
            name='logsig2_w',
            shape=(self.output_dim, self.input_dim),
            initializer='glorot_normal',
            trainable=True,
            dtype=K.floatx(),
        )
        if self.bayesian_bias:
            self.bias_logsig2 = self.add_weight(
                name='bias_logsig2',
                shape=(self.output_dim,),
                initializer='glorot_normal',
                trainable=True,
                dtype=K.floatx(),
            )

        # Draw the frozen eval samples once at build time, so the replica is ready to eval immediately after building.
        self.random = self.add_weight(
            name='random_w',
            shape=(self.output_dim, self.input_dim),
            initializer='zeros',
            trainable=False,
            dtype=K.floatx(),
        )
        if self.bayesian_bias:
            self.random_b = self.add_weight(
                name='random_b',
                shape=(self.output_dim,),
                initializer='zeros',
                trainable=False,
                dtype=K.floatx(),
            )

        self.reset_parameters()
        self.reset_random()

    def reset_parameters(self):
        stdv = 1.0 / math.sqrt(self.input_dim)
        self.bias.assign(Kops.zeros_like(self.bias))
        self.mu_w.assign(krandom.normal(self.mu_w.shape, mean=0.0, stddev=stdv, dtype=K.floatx()))
        self.logsig2_w.assign(
            krandom.normal(self.logsig2_w.shape, mean=self.std_init, stddev=0.001, dtype=K.floatx())
        )
        if self.bayesian_bias:
            # Start the bias posterior near-deterministic, like the weights.
            self.bias_logsig2.assign(
                krandom.normal(
                    self.bias_logsig2.shape, mean=self.std_init, stddev=0.001, dtype=K.floatx()
                )
            )

    @property
    def s2_w(self):
        """Weight variance, sigma_w^2 = exp(logsig2_w)."""
        return Kops.exp(Kops.clip(self.logsig2_w, self.lbound, self.ubound))

    @property
    def s2_b(self):
        """Bias variance, sigma_b^2 = exp(bias_logsig2)."""
        return Kops.exp(Kops.clip(self.bias_logsig2, self.lbound, self.ubound))

    def enable_map(self):
        self.map = True

    def disable_map(self):
        self.map = False

    def reset_random(self):
        """Redraw the frozen eval samples. Optional - build() already draws once;
        call this only if you want a fresh posterior sample for the replica."""
        self.random.assign(krandom.normal(self.random.shape, dtype=K.floatx()))
        if self.bayesian_bias:
            self.random_b.assign(krandom.normal(self.random_b.shape, dtype=K.floatx()))
        self.map = False

    def train(self):
        self.training = True

    def eval(self):
        self.training = False

    def _gaussian_kl(self, mu, logsig2):
        """Analytic KL[ N(mu, e^logsig2) || N(0, 1/prior_prec) ], summed."""
        logsig2 = Kops.clip(logsig2, self.lbound, self.ubound)
        return 0.5 * Kops.sum(
            self.prior_prec * (Kops.square(mu) + Kops.exp(logsig2))
            - logsig2
            - 1.0
            - Kops.log(self.prior_prec)
        )

    def kl_loss(self):
        kl = self._gaussian_kl(self.mu_w, self.logsig2_w)
        if self.bayesian_bias:
            kl += self._gaussian_kl(self.bias, self.bias_logsig2)
        return kl

    def _forward_map_inference(self, input):
        """Deterministic: posterior means for both weights and bias."""
        # use maximum-a-posteriori (MAP) method for computing predictions
        # sets each weight to its mean value (turns off weight sampling)
        return Kops.matmul(input, Kops.transpose(self.mu_w)) + self.bias

    def _forward_sample_activations(self, input):
        """
        Local reparameterization trick (training). Sample the pre-activation
        directly from its implied Gaussian. The bias variance (if Bayesian) is
        added straight into the activation variance. LRT is more efficient
        and leads to an estimate of the gradient with smaller variance.
        https://arxiv.org/pdf/1506.02557.pdf
        """
        act_mu = Kops.matmul(input, Kops.transpose(self.mu_w)) + self.bias
        act_var = Kops.matmul(Kops.square(input), Kops.transpose(self.s2_w))
        if self.bayesian_bias:
            act_var = act_var + self.s2_b
        act_var = act_var + self.eps
        noise = krandom.normal(Kops.shape(act_mu), dtype=input.dtype)
        return act_mu + Kops.sqrt(act_var) * noise

    def _forward_sample_weights(self, input):
        """
        Standard reparameterization (eval). Draw weights (and, if Bayesian, the
        bias) once, caching the standard normals as plain attributes so the
        replica stays frozen until reset_random().
        """
        weight = self.mu_w + Kops.sqrt(self.s2_w) * self.random
        bias = self.bias
        if self.bayesian_bias:
            bias = self.bias + Kops.sqrt(self.s2_b) * self.random_b
        return Kops.matmul(input, Kops.transpose(weight)) + bias

    def call(self, input):
        if self.training:
            return self._forward_sample_activations(input)
        if self.map:
            return self._forward_map_inference(input)
        return self._forward_sample_weights(input)


class Dense(KerasDense, MetaLayer):
    def __init__(self, **kwargs):
        # Set default dtype to tf.float64 if not provided
        if 'dtype' not in kwargs:
            kwargs['dtype'] = tf.float64
        super().__init__(**kwargs)


def dense_per_flavour(basis_size=8, kernel_initializer="glorot_normal", **dense_kwargs):
    """
    Generates a list of layers which can take as an input either one single layer
    or a list of the same size
    If taking one single layer, this one single layer will be the input of every layer in the list.
    If taking a list of layer of the same size, each layer on the list will take
    as input the layer on the input list in the same position.

    Note that, if the initializer is seeded, it should be a list where the seed is different
    for each element.

    i.e., if `basis_size` is 3 and is taking as input one layer A the output will be:
        [B1(A), B2(A), B3(A)]
    if taking, instead, a list [A1, A2, A3] the output will be:
        [B1(A1), B2(A2), B3(A3)]
    """
    if isinstance(kernel_initializer, str):
        kernel_initializer = basis_size * [kernel_initializer]

    # Need to generate a list of dense layers
    dense_basis = [
        base_layer_selector("dense", kernel_initializer=initializer, **dense_kwargs)
        for initializer in kernel_initializer
    ]

    def apply_dense(xinput):
        """
        The input can be either one single layer or a list of layers of
        length `basis_size`

        If taking one single layer, this one single layer will be the input of every
        layer in the list.
        If taking a list of layer of the same size, each layer on the list will take
        as input the layer on the input list in the same position.
        """
        if isinstance(xinput, (list, tuple)):
            if len(xinput) != basis_size:
                raise ValueError(f"""The input of the dense_per_flavour and the basis_size
doesn't match, got a list of length {len(xinput)} for a basis_size of {basis_size}""")
            results = [dens(ilayer) for dens, ilayer in zip(dense_basis, xinput)]
        else:
            results = [dens(xinput) for dens in dense_basis]

        return results

    return apply_dense


layers = {
    "dense": (
        Dense,
        {
            "kernel_initializer": "glorot_normal",
            "units": 5,
            "activation": "sigmoid",
            "kernel_regularizer": None,
            "dtype": tf.float64,
        },
    ),
    "dense_per_flavour": (
        dense_per_flavour,
        {
            "kernel_initializer": "glorot_normal",
            "units": 5,
            "activation": "sigmoid",
            "basis_size": 8,
            "dtype": tf.float64,
        },
    ),
    "LSTM": (
        LSTM_modified,
        {"kernel_initializer": "glorot_normal", "units": 5, "activation": "sigmoid"},
    ),
    "VBDense": (
        VBDense,
        {
            "in_features": None,
            "out_features": None,
            "prior_prec": None,
            "std_init": None,
            "bayesian_bias": False,
            "map": False,
        },
    ),
    "dropout": (Dropout, {"rate": 0.0}),
    "concatenate": (Concatenate, {}),
}

regularizers = {'l1_l2': (l1_l2, {'l1': 0.0, 'l2': 0.0})}


def base_layer_selector(layer_name, **kwargs):
    """
    Given a layer name, looks for it in the `layers` dictionary and returns an instance.

    The layer dictionary defines a number of defaults
    but they can be overwritten/enhanced through kwargs

    Parameters
    ----------
        `layer_name`
            str with the name of the layer
        `**kwargs`
            extra optional arguments to pass to the layer (beyond their defaults)
    """
    try:
        layer_tuple = layers[layer_name]
    except KeyError as e:
        raise NotImplementedError(
            f"Layer not implemented in keras_backend/base_layers.py: {layer_name}"
        ) from e

    layer_class = layer_tuple[0]
    layer_args = layer_tuple[1]

    for key, value in kwargs.items():
        # Check whether the activation function is a custom one
        if key == "activation":
            value = custom_activations.get(value, value)
        if key in layer_args.keys():
            layer_args[key] = value
        if key == "name":
            layer_args[key] = value

    return layer_class(**layer_args)


def regularizer_selector(reg_name, **kwargs):
    """Given a regularizer name looks in the `regularizer` dictionary and
    return an instance.

    The regularizer dictionary defines defaults for regularizers but these can
    be overwritten by supplying kwargs

    Parameters
    ----------
    layer_name
        str with the name of the regularizer
    **kwargs
        extra optional arguments to pass to the regularizer

    """
    if reg_name is None:
        return None

    try:
        reg_tuple = regularizers[reg_name]
    except KeyError:
        raise NotImplementedError(
            f"Regularizer not implemented in keras_backend/base_layers.py: {reg_name}"
        )

    reg_class = reg_tuple[0]
    reg_args = reg_tuple[1]

    for key, value in kwargs.items():
        if key in reg_args.keys():
            reg_args[key] = value

    return reg_class(**reg_args)
