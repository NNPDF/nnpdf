"""
    Module containg the losses to be apply to the models as layers

    The layer take the input from the model and acts on it producing a score function.
    For instance, in the case of the chi2 (``LossInvcovmat``) the function takes only
    the prediction of the model and, during instantiation, took the real data to compare with
    and the covmat.

"""

from unittest import result
import numpy as np
from keras import backend as K
import logging

from keras import ops as Kops
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

    def build(self, input_shape):
        """Transform the inverse covmat and the mask into
        weights of the layers"""
        init = MetaLayer.init_constant(self._invcovmat)
        self.kernel = self.builder_helper("invcovmat", self._invcovmat.shape, init, trainable=False)
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
        obs_diff_raw = self._y_true - y_pred
        # TODO: most of the time this is a y * I multiplication and can be skipped
        # benchmark how much time (if any) is lost in this in actual fits for the benefit of faster kfolds
        obs_diff = op.op_multiply([obs_diff_raw, self.mask])

        # The experimental loss doesn't depend on replicas, so it doesn't have a replica axis and
        # must be treated separately
        experimental_loss = len(self.kernel.shape) == 2
        one_replica = obs_diff.shape[1] == 1

        if one_replica:  # einsum is not well suited for CPU, so use tensordot if single replica
            kernel = self.kernel if experimental_loss else self.kernel[0]
            right_dot = op.tensor_product(kernel, obs_diff[0, 0, :], axes=1)
            loss = op.tensor_product(obs_diff[0, :, :], right_dot, axes=1)
        else:
            einstr = "bri, ij, brj -> r" if experimental_loss else "bri, rij, brj -> r"
            loss = op.einsum(einstr, obs_diff, self.kernel, obs_diff)

        return loss
    
class LossKL(MetaLayer):
    """Adds the KL divergence term once to the total training loss."""

    def __init__(self, vb_layers, kl_beta, **kwargs):
        self.vb_layers = vb_layers
        self.kl_beta = kl_beta
        super().__init__(**kwargs)

    def call(self, y_pred, **kwargs):
        kl = K.constant(0.0, dtype=K.floatx())
        for layer in self.vb_layers:
            kl += Kops.cast(layer.kl_loss(), K.floatx())
        result = kl * Kops.cast(self.kl_beta.read_value(), K.floatx())
        return op.reshape(result, (1,))


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
        loss = op.elu(-y_pred, alpha=self.alpha)
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
        y = y_pred * y_pred
        # Sum over the batch and the datapoints
        return op.sum(y, axis=[0, -1])
    

class LossRepulsion(MetaLayer):
    r"""
    Function-space SVGD repulsion across replicas (Repulsive Deep Ensembles).

    This is the n3fit translation of the per-batch repulsive-ensemble term:

      * the kernel is built in *function space*, i.e. directly on the PDF values
        evaluated on an anchor x-grid (the training grid), rather than on the
        per-data-point loss profile.  This is the prescription of
        "Repulsive Deep Ensembles are Bayesian" (D'Angelo & Fortuin 2021) and is
        the natural object in n3fit because the full multi-replica PDF tensor is
        already available in the graph.

      * the repulsion is NOT divided by the number of data points.  The
        prefactor is just ``beta``.

    Given the full multi-replica PDF on the anchor grid, shape ``(1, R, Nx, Fl)``,
    we flatten each replica to a vector ``f_r in R^{Nx*Fl}`` and build the RBF
    kernel:

        K_{kj} = exp( - ||f_k - f_j||^2 / (2 * sigma) )
        sigma  = median_{kj} ||f_k - f_j||^2 / (2 * log(R + 1))      (Liu & Wang 2016)

    The term added to the (already replica-summed) training loss is the
    gradient-only SVGD:

        L_rep = beta * sum_k ( rowsum_k / stop_grad(rowsum_k) - 1 ),
        rowsum_k = sum_j K_{kj}

    Numerically ``L_rep == 0`` (each summand is 1 - 1), so the *reported* chi2,
    the validation chi2 and therefore early stopping are left untouched, while
    the gradient is exactly the repulsive force that pushes replicas apart in
    function space.

    Parameters
    ----------
    beta : float
        Repulsion strength (the ``beta`` prefactor). 0 recovers a plain
        deep ensemble trained jointly.
    bandwidth : float or None
        If None (default) use the median heuristic each step (recomputed,
        detached). If a float is given it is used as a *fixed* ``sqrt(sigma)``
        (i.e. ``sigma = bandwidth**2``), which is sometimes preferable for very
        large R after a warm-up.
    flavour_weights : list/np.ndarray or None
        Optional per-flavour weights of length Fl (14 in the FK basis) applied to
        the function-space distance, e.g. to focus repulsion on the gluon and the
        quark singlet. If None, all flavours contribute equally.
    """

    def __init__(self, beta=1.0, bandwidth=None, flavour_weights=None, **kwargs):
        self.beta = float(beta)
        self.bandwidth = bandwidth
        if flavour_weights is not None:
            fw = np.array(flavour_weights, dtype=np.float64).reshape(1, 1, -1)
            self._flavour_weights = op.numpy_to_tensor(fw)
        else:
            self._flavour_weights = None
        super().__init__(**kwargs)

    def call(self, pdf_anchors, **kwargs):
        # pdf_anchors: (1, R, Nx, Fl)  -> drop the batch axis -> (R, Nx, Fl)
        x = pdf_anchors[0]
        if self._flavour_weights is not None:
            x = x * self._flavour_weights

        r = Kops.shape(x)[0]
        flat = op.reshape(x, (r, -1))               # (R, D), D = Nx*Fl
        flat_detached = Kops.stop_gradient(flat)    # gradient only via `flat`

        # Pairwise squared distances ||f_k - f_j||^2, gradient through f_k only.
        # diff[k, j, :] = f_k - stop_grad(f_j)
        diff = Kops.expand_dims(flat, 1) - Kops.expand_dims(flat_detached, 0)  # (R, R, D)
        dnorm2 = op.sum(Kops.square(diff), axis=2)                             # (R, R)

        # Bandwidth (median estimator), detached so it is treated as a constant.
        if self.bandwidth is None:
            med = Kops.quantile(Kops.stop_gradient(dnorm2), 0.5)
            denom = 2.0 * op.op_log(Kops.cast(r, dnorm2.dtype) + 1.0)
            sigma = med / denom
            sigma = Kops.maximum(sigma, Kops.cast(1e-12, dnorm2.dtype))
        else:
            sigma = Kops.cast(self.bandwidth ** 2, dnorm2.dtype)

        kmat = Kops.exp(-dnorm2 / (2.0 * sigma))       # (R, R)
        rowsum = op.sum(kmat, axis=1)                  # (R,)

        # SVGD gradient-only trick: value == 0, gradient == repulsive force.
        rep = op.sum(rowsum / Kops.stop_gradient(rowsum) - 1.0)
        result = Kops.cast(self.beta, rep.dtype) * rep
        return op.reshape(result, (1,))