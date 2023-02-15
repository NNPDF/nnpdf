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

    Both the inverse covmat and the mask (if any) are stored as layer weights
    and can be updated at any points either directly or by using the
    ``update_mask`` and ``add_covmat`` methods.

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

    def __init__(self, invcovmat, y_true, mask=None, covmat=None, **kwargs):
        # If we have a diagonal matrix, padd with 0s and hope it's not too heavy on memory
        if len(invcovmat.shape) == 1:
            invcovmat = np.diag(invcovmat)
        self._invcovmat = op.numpy_to_tensor(invcovmat)
        self._covmat = covmat
        self._y_true = op.numpy_to_tensor(y_true)
        self._ndata = y_true.shape[-1]
        if mask is None or all(mask):
            self._mask = None
        else:
            mask = np.array(mask, dtype=np.float32).reshape((1, 1, -1))
            self._mask = op.numpy_to_tensor(mask)
        super().__init__(**kwargs)

    def build(self):
        """Transform the inverse covmat and the mask into
        weights of the layers"""
        init = MetaLayer.init_constant(self._invcovmat)
        self.kernel = self.builder_helper(
            "invcovmat", self._invcovmat.shape, init, trainable=False
        )
        mask_shape = (1, 1, self._ndata)
        if self._mask is None:
            init_mask = MetaLayer.init_constant(np.ones(mask_shape))
        else:
            init_mask = MetaLayer.init_constant(self._mask)
        self.mask = self.builder_helper("mask", mask_shape, init_mask, trainable=False)

    def add_covmat(self, covmat):
        """Add a piece to the inverse covmat weights
        Note, however, that the _covmat attribute of the layer will
        still refer to the original data covmat
        """
        new_covmat = np.linalg.inv(self._covmat + covmat)
        self.kernel.assign(new_covmat)

    def update_mask(self, new_mask):
        """Update the mask"""
        self.mask.assign(new_mask)

    def call(self, y_pred, **kwargs):
        tmp_raw = self._y_true - y_pred
        # TODO: most of the time this is a y * I multiplication and can be skipped
        # benchmark how much time (if any) is lost in this in actual fits for the benefit of faster kfolds
        tmp = op.op_multiply([tmp_raw, self.mask])
        if tmp.shape[1] == 1:
            # einsum is not well suited for CPU, so use tensordot if not multimodel
            right_dot = op.tensor_product(self.kernel, tmp[0, 0, :], axes=1)
            res = op.tensor_product(tmp[0, :, :], right_dot, axes=1)
        else:
            if len(self.kernel.shape) == 3:
                res = op.einsum("bri, rij, brj -> r", tmp, self.kernel, tmp)
            else:
                res = op.einsum("bri, ij, brj -> r", tmp, self.kernel, tmp)
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
