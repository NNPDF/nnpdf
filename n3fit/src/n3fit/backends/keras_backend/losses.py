"""
    Module containing a list of loss functions availables to the fitting code
"""

import tensorflow as tf
from tensorflow.keras import backend as K


def l_invcovmat(invcovmat_np):
    """
    Returns a loss function such that:
    L = \sum_{ij} (yt - yp)_{i} invcovmat_{ij} (yt - yp)_{j}
    """
    invcovmat = K.constant(invcovmat_np)

    def true_loss(y_true, y_pred):
        """
            y_true: (1, N)
            y_pred: (1, replicas, N)
        """
        # (yt - yp) * covmat * (yt - yp)
        tmp = y_true - y_pred
        return tf.einsum('bri,ij,brj->b', tmp, invcovmat, tmp)

    return true_loss


def l_positivity(alpha=1e-7):
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

    def true_loss(y_true, y_pred):
        y = -y_pred
        loss = K.elu(y, alpha=alpha)
        res = K.sum(loss)
        return tf.reshape(res, (-1,))

    return true_loss


def l_integrability():
    """
    Returns (y_pred)*(y_pred)
    """

    def true_loss(y_true, y_pred):
        loss = K.square(y_pred)
        res = K.sum(loss, keepdims=True)
        return tf.reshape(res, (-1,))

    return true_loss


def l_diaginvcovmat(diaginvcovmat_np):
    """
    Returns a loss function such that:
    L = sum_{i} (yt - yp)_{i} invcovmat_{ii} (yt - yp)_{i}
    diaginvcovmat_np should be 1d
    """
    invcovmat = K.constant(diaginvcovmat_np)

    def true_loss(y_true, y_pred):
        tmp = y_true - y_pred
        res = tf.tensordot(invcovmat, K.transpose(tmp * tmp), axes=1)
        return tf.reshape(res, (-1,))

    return true_loss
