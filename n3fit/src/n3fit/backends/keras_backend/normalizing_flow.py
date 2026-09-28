"""
RealNVP normalizing flow used to enrich VBDense's mean-field Gaussian posterior.

An identity-initialised stack of affine coupling layers that maps a 
mean-field Gaussian sample ``z`` into a (possibly multimodal)
weight sample ``w = flow(z)``, with a tractable per-sample log-det Jacobian so that
``log q(w) = log N(z; mu, sigma^2) - log|det d(flow)/dz|`` stays exact and can be used
as a MC entropy term in a KL estimate (see ``VBDense.kl_loss`` in
``base_layers.py``).
"""

import numpy as np
import keras.backend as K
from keras import ops as Kops
from keras.initializers import Constant
from keras.layers import Dense, Layer


class RealNVPCouplingLayer(Layer):
    """
    Affine coupling: with binary mask m (conditioning half) and im = 1 - m,
        y = m*x + im*(x*exp(s(x*m)) + t(x*m))
        log|det J| = sum_dim( im * s )          # per sample, shape (batch,)
    s2, t2 zero-initialised -> the layer starts as the identity.
    """

    def __init__(self, mask, hidden_units=24, scale_clip=2.0, **kwargs):
        super().__init__(**kwargs)
        self._mask_np = np.asarray(mask, dtype=K.floatx())
        self.hidden_units = hidden_units
        self.scale_clip = scale_clip

    def build(self, input_shape):
        d = int(input_shape[-1])
        dtype = K.floatx()
        self.s1 = Dense(
            self.hidden_units,
            activation='relu',
            kernel_initializer='glorot_normal',
            name=f'{self.name}_s1',
            dtype=dtype,
        )
        self.s2 = Dense(
            d, activation='tanh', kernel_initializer='zeros', name=f'{self.name}_s2', dtype=dtype
        )
        self.t1 = Dense(
            self.hidden_units,
            activation='relu',
            kernel_initializer='glorot_normal',
            name=f'{self.name}_t1',
            dtype=dtype,
        )
        self.t2 = Dense(d, kernel_initializer='zeros', name=f'{self.name}_t2', dtype=dtype)
        self._mask = self.add_weight(
            name='mask',
            shape=(d,),
            initializer=Constant(self._mask_np),
            trainable=False,
            dtype=dtype,
        )
        super().build(input_shape)

    def call(self, x):
        m = self._mask
        im = 1.0 - m
        s = self.scale_clip * self.s2(self.s1(x * m))
        t = self.t2(self.t1(x * m))
        y = m * x + im * (x * Kops.exp(s) + t)
        log_det = Kops.sum(im * s, axis=-1)  # (batch,)
        return y, log_det


class RealNVPFlow(Layer):
    """Stack of alternating-mask RealNVP couplings. call(x) -> (y, log_det[per-sample])."""

    def __init__(self, n_flows, dim, hidden_units=24, scale_clip=2.0, **kwargs):
        super().__init__(**kwargs)
        if dim < 2:
            raise ValueError(f"RealNVPFlow: dim must be >= 2, got {dim}")
        masks = []
        for i in range(n_flows):
            m = np.zeros(dim, dtype=K.floatx())
            m[: dim // 2] = 1.0 if i % 2 == 0 else 0.0
            m[dim // 2 :] = 0.0 if i % 2 == 0 else 1.0
            masks.append(m)
        self.coupling_layers = [
            RealNVPCouplingLayer(
                m, hidden_units=hidden_units, scale_clip=scale_clip, name=f'{self.name}_c{i}'
            )
            for i, m in enumerate(masks)
        ]

    def call(self, x):
        log_det = Kops.zeros(Kops.shape(x)[0], dtype=x.dtype)
        y = x
        for cl in self.coupling_layers:
            y, ld = cl(y)
            log_det = log_det + ld
        return y, log_det
