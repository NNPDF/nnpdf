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
from .normalizing_flow import RealNVPFlow
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
        std_init: float = None,
        map: bool = False,
        bayesian_bias: bool = False,
        use_flow: bool = False,
        n_flows: int = 3,
        flow_hidden: int = 24,
    ):
        super().__init__()
        self.output_dim = out_features
        self.input_dim = in_features
        self.map = map
        self.prior_prec = float(prior_prec)
        if std_init is None:
            self.std_init = -math.log(self.prior_prec)
        else:
            self.std_init = float(std_init)
        self.bayesian_bias = bayesian_bias
        self.lbound = -30 if K.floatx() == 'float64' else -20
        self.ubound = 11
        self.eps = 1e-12 if K.floatx() == 'float64' else 1e-8
        self.training = True
        # Frozen eval draws (plain attributes, NOT tracked weights).
        self.random = None
        self.random_b = None
        # Optional RealNVP flow on top of the mean-field Gaussian base (see
        # normalizing_flow.py). Off by default: plain VBDense with use_flow=False
        self.use_flow = bool(use_flow)
        self.n_flows = int(n_flows)
        self.flow_hidden = int(flow_hidden)
        self.flow = None
        # z/w/logdet cached by the most recent stochastic forward pass, read back by
        # kl_loss() for the flow's Monte-Carlo KL estimate (see model_trainer.py's
        # _pdf_injection for why this ordering is load-bearing)
        self._last_z = None
        self._last_w = None
        self._last_logdet = None

    def build(self, input_shape):
        self.bias = self.add_weight(
            name='bias',
            shape=(self.output_dim,),
            initializer='glorot_normal',
            trainable=True,
            dtype=K.floatx(),
        )
        weight_shape = (self.output_dim, self.input_dim)
        self.mu_w = self.add_weight(
            name='mu_w',
            shape=weight_shape,
            initializer='glorot_normal',
            trainable=True,
            dtype=K.floatx(),
        )
        self.logsig2_w = self.add_weight(
            name='logsig2_w',
            shape=weight_shape,
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

        # Training noise: one standard-normal draw per training step, shared by all calls of
        # the layer within the step (all FK x-grids, x=1, sum-rule grid), so that every call
        # sees the same weights. Redrawn by callbacks.ResampleTrainNoise. Plain tf.Variables:
        # not trainable and not part of layer.weights, so the saved weights are unchanged.
        self.train_noise = tf.Variable(
            tf.zeros(weight_shape, dtype=K.floatx()), trainable=False, name='train_noise_w'
        )
        self.train_noise_eps = None
        if self.bayesian_bias:
            self.train_noise_b = tf.Variable(
                tf.zeros((self.output_dim,), dtype=K.floatx()), trainable=False, name='train_noise_b'
            )

        if self.use_flow:
            dim = self.output_dim * self.input_dim
            self.flow = RealNVPFlow(
                self.n_flows, dim, hidden_units=self.flow_hidden, name=f'{self.name}_flow'
            )
            # Warm up so the flow's variables are tracked from the start (identity at init)
            _ = self.flow(Kops.zeros((1, dim), dtype=K.floatx()))

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
        call this only if you want a fresh posterior sample for the replica.
        """
        self.random.assign(krandom.normal(self.random.shape, dtype=K.floatx()))
        if self.bayesian_bias:
            self.random_b.assign(krandom.normal(self.random_b.shape, dtype=K.floatx()))
        self.map = False

    def resample_train_noise(self, rng):
        """New training-noise draw from the numpy Generator rng (once per training step)."""
        for var in self.train_noise_variables():
            var.assign(rng.standard_normal(var.shape).astype(var.dtype.as_numpy_dtype))

    def train_noise_variables(self):
        out = [self.train_noise]
        if self.train_noise_eps is not None:
            out.append(self.train_noise_eps)
        if self.bayesian_bias:
            out.append(self.train_noise_b)
        return out

    def train(self):
        self.training = True

    def eval(self):
        self.training = False

    def _gaussian_kl(self, mu, logsig2):
        """Analytic KL[ N(mu, e^logsig2) || N(0, 1/prior_prec) ], summed."""
        logsig2 = Kops.clip(logsig2, self.lbound, self.ubound)
        # NOTE: Kops.log(python_float) ignores keras.config.floatx() and returns a
        # float32 tensor regardless (a Keras-3 keras.ops quirk), which crashes this
        # subtraction against a float64 logsig2 with a dtype mismatch. math.log is a
        # plain Python float and combines with the tensor via normal promotion.
        return 0.5 * Kops.sum(
            self.prior_prec * (Kops.square(mu) + Kops.exp(logsig2))
            - logsig2
            - 1.0
            - math.log(self.prior_prec)
        )

    def _log_prior(self, w):
        """log N(w; 0, 1/prior_prec), summed over all elements of w."""
        return -0.5 * Kops.sum(
            self.prior_prec * Kops.square(w) + math.log(2.0 * math.pi) - math.log(self.prior_prec)
        )

    def _mc_kl_flow(self):
        """
        Single-sample MC estimate of KL[q(w)||p(w)] = log q(w) - log p(w), for
        the flow-transformed posterior q(w) = flow(z), z ~ N(mu_w, sigma_w^2).

        Uses z/w/logdet cached by the most recent stochastic forward pass (there is no
        closed form once w is flow-transformed). At flow == identity this reduces, in
        expectation over z, to the analytic `_gaussian_kl(mu_w, logsig2_w)` above 
        """
        if self._last_z is None:
            raise RuntimeError(
                "VBDense.kl_loss() called with use_flow=True before any stochastic "
                "forward pass; the z/w/logdet cache is empty."
            )
        logsig2 = Kops.clip(self.logsig2_w, self.lbound, self.ubound)
        eps = (self._last_z - self.mu_w) / Kops.exp(0.5 * logsig2)
        log_q0 = -0.5 * Kops.sum(Kops.square(eps) + math.log(2.0 * math.pi) + logsig2)
        log_q = log_q0 - Kops.sum(self._last_logdet)
        log_p = self._log_prior(self._last_w)
        return log_q - log_p

    def kl_loss(self):
        if self.use_flow:
            kl = self._mc_kl_flow()
        else:
            kl = self._gaussian_kl(self.mu_w, self.logsig2_w)
        if self.bayesian_bias:
            kl += self._gaussian_kl(self.bias, self.bias_logsig2)
        return kl

    def _sample_weight_and_logq(self, z):
        """
        Push a mean-field sample z (shape (out,in)) through the flow, if attached.
        Returns (w, logdet) with w shaped like z and logdet a per-sample tensor
        (shape (1,), since VBDense always draws exactly one weight sample at a time).
        With use_flow=False this is the identity: w=z, logdet=0.
        """
        if not self.use_flow:
            return z, Kops.zeros((1,), dtype=z.dtype)
        dim = self.output_dim * self.input_dim
        z_flat = Kops.reshape(z, (1, dim))
        w_flat, logdet = self.flow(z_flat)
        w = Kops.reshape(w_flat, (self.output_dim, self.input_dim))
        return w, logdet

    def _forward_map_inference(self, input):
        """Deterministic: posterior means for both weights and bias."""
        # use maximum-a-posteriori (MAP) method for computing predictions
        # sets each weight to its mean value (turns off weight sampling)
        weight = self.mu_w
        return Kops.matmul(input, Kops.transpose(weight)) + self.bias

    def _forward_sample_weights_train(self, input):
        """
        Training path. One weight sample per training step, w = mu + sigma * train_noise,
        shared by all x points and all calls of the layer within the step. The local
        reparameterization trick is not used: it draws independent noise per x point,
        whereas the FK convolutions and the chi2 couple all x points, so the expected chi2
        depends on the covariance of the PDF between different x points.
        """
        weight = self.mu_w + Kops.sqrt(self.s2_w) * self.train_noise
        bias = self.bias
        if self.bayesian_bias:
            bias = self.bias + Kops.sqrt(self.s2_b) * self.train_noise_b
        return Kops.matmul(input, Kops.transpose(weight)) + bias


    def _forward_sample_weights(self, input):
        """
        Standard reparameterization (eval). Draw weights (and, if Bayesian, the
        bias) once, caching the standard normals as plain attributes so the
        replica stays frozen until reset_random(). If a flow is attached, the frozen
        Gaussian sample is pushed through it too, so eval/pseudo-replica sampling uses
        the same trained-Gaussian-plus-flow posterior as training.
        """
        z = self.mu_w + Kops.sqrt(self.s2_w) * self.random
        weight, _ = self._sample_weight_and_logq(z)
        bias = self.bias
        if self.bayesian_bias:
            bias = self.bias + Kops.sqrt(self.s2_b) * self.random_b
        return Kops.matmul(input, Kops.transpose(weight)) + bias

    def _forward_sample_weights_stochastic(self, input):
        """
        Explicit reparameterization (training, flow-enabled layers only).

        LRT (`_forward_sample_activations`) cannot be used here: it propagates the
        pre-activation distribution analytically, which requires q(w) to be Gaussian --
        no longer true once w = flow(z). Instead, draw a FRESH weight sample every call
        (unlike the frozen `self.random` used at eval) so each training step is an
        independent MC draw of the reverse-KL objective, and cache z/w/logdet
        for kl_loss() to reuse (see _mc_kl_flow).
        """
        z = self.mu_w + Kops.sqrt(self.s2_w) * self.train_noise
        w, logdet = self._sample_weight_and_logq(z)

        bias = self.bias
        if self.bayesian_bias:
            bias = self.bias + Kops.sqrt(self.s2_b) * self.train_noise_b

        self._last_z, self._last_w, self._last_logdet = z, w, logdet

        return Kops.matmul(input, Kops.transpose(w)) + bias

    def call(self, input):
        if self.training:
            if self.use_flow:
                return self._forward_sample_weights_stochastic(input)
            return self._forward_sample_weights_train(input)
        if self.map:
            return self._forward_map_inference(input)
        return self._forward_sample_weights(input)

class CorrelatedLowRankVBDense(Layer):
    """
    Variational Bayesian dense layer with a correlated weight posterior.

        q(w) = N(mu, Sigma),   Sigma = D + U U^T,
        D = diag(exp(logsig2)) > 0,   U in R^{N x rank},

    where N = out_features * in_features (+ out_features if the bias is
    Bayesian). Sigma is positive definite and invertible by construction,
    costs O(N * rank) parameters instead of O(N^2), and its log-determinant
    follows from the rank x rank capacitance matrix (matrix determinant
    lemma). Setting rank=0 reproduces `VBDense` exactly.

    Backend-agnostic: written with `keras.ops` / `keras.random`.
    """

    def __init__(
        self,
        out_features: int,
        in_features: int,
        rank: int = 4,
        prior_prec: float = 1.0,
        std_init: float = None,
        u_init: float = 1e-4,
        map: bool = False,
        bayesian_bias: bool = False,
    ):
        super().__init__()
        if rank < 0:
            raise ValueError("rank must be non-negative")
        self.output_dim = out_features
        self.input_dim = in_features
        self.rank = int(rank)
        self.map = map
        self.prior_prec = float(prior_prec)
        if std_init is None:
            self.std_init = -math.log(self.prior_prec)
        else:
            self.std_init = float(std_init)
        self.u_init = float(u_init)
        self.bayesian_bias = bayesian_bias
        self.lbound = -30 if K.floatx() == 'float64' else -20
        self.ubound = 11
        self.eps = 1e-12 if K.floatx() == 'float64' else 1e-8
        self.training = True
        # Frozen eval draws, assigned in build().
        self.random = None
        self.random_b = None
        self.random_eps = None

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
        if self.rank > 0:
            self.u_w = self.add_weight(
                name='u_w',
                shape=(self.rank, self.output_dim, self.input_dim),
                initializer='zeros',
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
            if self.rank > 0:
                self.u_b = self.add_weight(
                    name='u_b',
                    shape=(self.rank, self.output_dim),
                    initializer='zeros',
                    trainable=True,
                    dtype=K.floatx(),
                )

        # Frozen eval draws (non-trainable), so the replica is ready to eval
        # immediately after building.
        self.random = self.add_weight(
            name='random_w',
            shape=(self.output_dim, self.input_dim),
            initializer='zeros',
            trainable=False,
            dtype=K.floatx(),
        )
        if self.rank > 0:
            self.random_eps = self.add_weight(
                name='random_eps',
                shape=(self.rank,),
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

        # Training noise, one draw per training step (see VBDense.build)
        self.train_noise = tf.Variable(
            tf.zeros((self.output_dim, self.input_dim), dtype=K.floatx()),
            trainable=False,
            name='train_noise_w',
        )
        self.train_noise_eps = None
        if self.rank > 0:
            self.train_noise_eps = tf.Variable(
                tf.zeros((self.rank,), dtype=K.floatx()), trainable=False, name='train_noise_eps'
            )
        if self.bayesian_bias:
            self.train_noise_b = tf.Variable(
                tf.zeros((self.output_dim,), dtype=K.floatx()), trainable=False, name='train_noise_b'
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
        if self.rank > 0:
            # Start essentially mean-field; the ELBO has to grow the
            # correlated directions.
            self.u_w.assign(
                krandom.normal(self.u_w.shape, mean=0.0, stddev=self.u_init, dtype=K.floatx())
            )
        if self.bayesian_bias:
            self.bias_logsig2.assign(
                krandom.normal(
                    self.bias_logsig2.shape, mean=self.std_init, stddev=0.001, dtype=K.floatx()
                )
            )
            if self.rank > 0:
                self.u_b.assign(
                    krandom.normal(self.u_b.shape, mean=0.0, stddev=self.u_init, dtype=K.floatx())
                )

    @property
    def s2_w(self):
        """Diagonal weight variance, exp(logsig2_w)."""
        return Kops.exp(Kops.clip(self.logsig2_w, self.lbound, self.ubound))

    @property
    def s2_b(self):
        """Diagonal bias variance, exp(bias_logsig2)."""
        return Kops.exp(Kops.clip(self.bias_logsig2, self.lbound, self.ubound))

    def enable_map(self):
        self.map = True

    def disable_map(self):
        self.map = False

    def reset_random(self):
        """Redraw the frozen eval sample (diagonal and low-rank parts)."""
        self.random.assign(krandom.normal(self.random.shape, dtype=K.floatx()))
        if self.rank > 0:
            self.random_eps.assign(krandom.normal(self.random_eps.shape, dtype=K.floatx()))
        if self.bayesian_bias:
            self.random_b.assign(krandom.normal(self.random_b.shape, dtype=K.floatx()))
        self.map = False

    def resample_train_noise(self, rng):
        """New training-noise draw from the numpy Generator rng (once per training step)."""
        for var in self.train_noise_variables():
            var.assign(rng.standard_normal(var.shape).astype(var.dtype.as_numpy_dtype))

    def train_noise_variables(self):
        out = [self.train_noise]
        if self.train_noise_eps is not None:
            out.append(self.train_noise_eps)
        if self.bayesian_bias:
            out.append(self.train_noise_b)
        return out

    def train(self):
        self.training = True

    def eval(self):
        self.training = False

    def _flat_posterior(self):
        """Flattened (mu, d, U) with layout [weights, (bias)]."""
        mu = Kops.reshape(self.mu_w, (-1,))
        d = Kops.reshape(self.s2_w, (-1,))
        u = None
        if self.rank > 0:
            u = Kops.transpose(Kops.reshape(self.u_w, (self.rank, -1)))
        if self.bayesian_bias:
            mu = Kops.concatenate([mu, Kops.reshape(self.bias, (-1,))], axis=0)
            d = Kops.concatenate([d, Kops.reshape(self.s2_b, (-1,))], axis=0)
            if self.rank > 0:
                u_b = Kops.transpose(Kops.reshape(self.u_b, (self.rank, -1)))
                u = Kops.concatenate([u, u_b], axis=0)
        return mu, d, u

    def kl_loss(self):
        """Analytic KL[ N(mu, D + U U^T) || N(0, 1/prior_prec) ]."""
        mu, d, u = self._flat_posterior()
        n = self.output_dim * self.input_dim
        if self.bayesian_bias:
            n += self.output_dim

        trace = Kops.sum(d)
        logdet = Kops.sum(Kops.log(d))

        if self.rank > 0:
            trace = trace + Kops.sum(Kops.square(u))
            # log|Sigma| = log|D| + log|I_R + U^T D^-1 U|
            g = u / Kops.expand_dims(Kops.sqrt(d), axis=-1)
            cap = Kops.eye(self.rank, dtype=K.floatx()) + Kops.matmul(Kops.transpose(g), g)
            chol = Kops.cholesky(cap)
            logdet = logdet + 2.0 * Kops.sum(Kops.log(Kops.diagonal(chol)))

        return 0.5 * (
            self.prior_prec * (trace + Kops.sum(Kops.square(mu)))
            - logdet
            - n
            - n * math.log(self.prior_prec)
        )

    def _forward_map_inference(self, input):
        """Deterministic: posterior means for both weights and bias."""
        return Kops.matmul(input, Kops.transpose(self.mu_w)) + self.bias

    def _forward_sample_weights_train(self, input):
        """
        Training path. One weight sample per training step, w = mu + D^(1/2) eta + U eps,
        with the standard normals of the current step (train_noise, train_noise_eps,
        train_noise_b), shared by all x points and all calls of the layer within the step.
        """
        weight = self.mu_w + Kops.sqrt(self.s2_w) * self.train_noise
        bias = self.bias
        if self.bayesian_bias:
            bias = bias + Kops.sqrt(self.s2_b) * self.train_noise_b
        if self.rank > 0:
            weight = weight + Kops.einsum("r,roi->oi", self.train_noise_eps, self.u_w)
            if self.bayesian_bias:
                bias = bias + Kops.einsum("r,ro->o", self.train_noise_eps, self.u_b)
        return Kops.matmul(input, Kops.transpose(weight)) + bias


    def _forward_sample_weights(self, input):
        """
        Eval path. One frozen weight draw, w = mu + D^(1/2) eta + U eps, with
        the standard normals held in non-trainable weights so the replica
        stays fixed until reset_random().
        """
        weight = self.mu_w + Kops.sqrt(self.s2_w) * self.random
        bias = self.bias
        if self.bayesian_bias:
            bias = bias + Kops.sqrt(self.s2_b) * self.random_b
        if self.rank > 0:
            weight = weight + Kops.einsum("r,roi->oi", self.random_eps, self.u_w)
            if self.bayesian_bias:
                bias = bias + Kops.einsum("r,ro->o", self.random_eps, self.u_b)
        return Kops.matmul(input, Kops.transpose(weight)) + bias

    def call(self, input):
        if self.training:
            return self._forward_sample_weights_train(input)
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
            "use_flow": False,
            "n_flows": 3,
            "flow_hidden": 24,
        },
    ),
    "VBDense_correlated": (
        CorrelatedLowRankVBDense,
        {
            "in_features": None,
            "out_features": None,
            "rank": 4,
            "prior_prec": None,
            "std_init": None,
            "u_init": 1e-4,
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
    layer_args = dict(layer_tuple[1])

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
