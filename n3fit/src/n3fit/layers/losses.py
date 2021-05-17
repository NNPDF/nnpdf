"""
    Module containg the losses to be apply to the models as layers

    The layer take the input from the model and acts on it producing a score function.
    For instance, in the case of the chi2 (``LossInvcovmat``) the function takes only
    the prediction of the model and, during instantiation, took the real data to compare with
    and the covmat.

"""
import numpy as np
from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


class LossInvcovmat(MetaLayer):
    """
    Loss function such that:
    L = \\sum_{ij} (yt - yp)_{i} invcovmat_{ij} (yt - yp)_{j}

    Takes as argument the inverse of the covmat and the target data.
    It also takes an optional argument to mask part of the predictions

    Example
    -------
    >>> import numpy as np
    >>> from n3fit.layers import losses
    >>> C = np.random.rand(5,5)
    >>> data = np.random.rand(1, 1, 5)
    >>> pred = np.random.rand(1, 1, 5)
    >>> invC = np.linalg.inv( C @ C.T)
    >>> loss_f = losses.LossInvcovmat(invC, data)
    >>> loss_f(pred).shape == 1
    True
    """

    def __init__(self, invcovmat, y_true, mask=None, **kwargs):
        # If we have a diagonal matrix, padd with 0s and hope it's not too heavy on memory
        if len(invcovmat.shape) == 1:
            invcovmat = np.diag(invcovmat)
        self.invcovmat = op.numpy_to_tensor(invcovmat)
        self.y_true = op.numpy_to_tensor(y_true)
        if mask is None or all(mask):
            self.mask = None
        else:
            mask = np.array(mask, dtype=np.float32).reshape((1,1,-1))
            self.mask = op.numpy_to_tensor(mask)
        super().__init__(**kwargs)

    def call(self, y_pred, **kwargs):
        tmp = self.y_true - y_pred
        if self.mask is not None:
            tmp = op.op_multiply([tmp, self.mask])
        if tmp.shape[1] == 1:
            # einsum is not well suited for CPU, so use tensordot if not multimodel
            right_dot = op.tensor_product(self.invcovmat, tmp[0,0,:], axes=1)
            res = op.tensor_product(tmp[0,:,:], right_dot, axes=1)
        else:
            res = op.einsum("bri, ij, brj -> r", tmp, self.invcovmat, tmp)
        return res


class LossLagrange(MetaLayer):
    """
    Abstract loss function to apply lagrange multipliers to a model.

        L = \\lambda * f(y)

    The form of f(y) is given by modifying the ``apply_loss`` method.
    It is possible to modify how the multiplication of the lambda factor is implemented
    by modifying the ``apply_multiplier`` method.

    The (non trainable) weight containing the multiplier is named ``lagMult``.
    """

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
    Returns L = \\lambda*elu(y_pred)

    The positivity loss is computed by inverting the sign of the
    datapoints and then applying the elu function, this function is
        f(x) = x if x > 0
        f(x) = alpha * (e^{x} - 1) if x < 0
    This is done to avoid a big discontinuity in the derivative at 0 when
    the lagrange multiplier is very big.
    In practice this function can produce results in the range (-alpha, inf)

    Example
    -------
    >>> import numpy as np
    >>> from n3fit.layers import losses
    >>> pred = np.random.rand(1, 1, 5)
    >>> alpha = 1e-7
    >>> c = 1e8
    >>> loss_f = losses.LossPositivity(c=c, alpha=alpha)
    >>> loss_f(pred) == -5*alpha
    True
    >>> loss_f(-pred) > c
    True
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

    Example
    -------
    >>> import numpy as np
    >>> from n3fit.layers import losses
    >>> pred = np.random.rand(1, 1, 5)
    >>> loss_f = losses.LossIntegrability(c=1e2)
    >>> loss_f(pred) > 0
    True
    """

    def apply_loss(self, y_pred):
        y = op.backend_function("square", y_pred)
        # Sum over the batch and the datapoints
        return op.sum(y, axis=[0, -1])
