"""
    Module containg the losses to be apply to the models as layers
"""

from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


class LossInvcovmat(MetaLayer):
    """
    Loss function such that:
    L = \sum_{ij} (yt - yp)_{i} invcovmat_{ij} (yt - yp)_{j}
    """

    def __init__(self, invcovmat, y_true, **kwargs):
        self.invcovmat = op.numpy_to_tensor(invcovmat)
        self.y_true = op.numpy_to_tensor(y_true)
        super().__init__(**kwargs)

    def call(self, y_pred, **kwargs):
        tmp = self.y_true - y_pred
        res = op.einsum("bri, ij, brj -> r", tmp, self.invcovmat, tmp)
        return res


class LossLagrange(MetaLayer):
    def __init__(self, c=1.0, **kwargs):
        self._initial_multiplier = c
        super().__init__(**kwargs)

    def build(self, input_shape):
        multiplier = MetaLayer.init_constant(self._initial_multiplier)
        self.kernel = self.builder_helper("lagMult", (1,), multiplier, trainable=False)
        super().build(input_shape)

    def apply_multiplier(self, y):
        return self.kernel * y

    def apply_loss(self, y):
        return y

    def call(self, y_pred, **kwargs):
        y = self.apply_multiplier(y_pred)
        return self.apply_loss(y)


class LossPositivity(LossLagrange):
    """
    Returns L = elu(y_pred) (considers y_true as 0)

    The positivity loss is computed by inverting the sign of the
    datapoints and then applying the elu function, this function is
        f(x) = x if x > 0
        f(x) = alpha * (e^{x} - 1) if x < 0
    This is done to avoid a big discontinuity in the derivative at 0 when
    the lagrange multiplier is very big.
    In practice this function can produce results in the range (-alpha, inf)
    """

    def __init__(self, alpha=1e-7, **kwargs):
        self.alpha = alpha
        super().__init__(**kwargs)

    def apply_loss(self, y_pred):
        loss = op.backend_function("elu", -y_pred, alpha=self.alpha)
        # Sum over the batch and the datapoints
        return op.sum(loss, axis=[0, -1])


class LossIntegrability(LossLagrange):
    """
    Returns L = (y_pred)*(y_pred)
    """

    def apply_loss(self, y_pred):
        y = op.backend_function("square", y_pred)
        # Sum over the batch and the datapoints
        return op.sum(y, axis=[0, -1])
