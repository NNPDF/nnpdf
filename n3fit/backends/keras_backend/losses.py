import keras.backend as K


def l_invcovmat(invcovmat_np):
    """
    Returns a loss function such that:
    L = \sum_{ij} (yt - yp)_{i} invcovmat_{ij} (yt - yp)_{j}
    """
    invcovmat = K.constant(invcovmat_np)

    def true_loss(y_true, y_pred):
        # (yt - yp) * covmat * (yt - yp)
        tmp = y_true - y_pred
        right_dot = K.tf.tensordot(invcovmat, K.transpose(tmp), axes=1)
        return K.tf.tensordot(tmp, right_dot, axes=1)

    return true_loss


def l_positivity(alpha=1e-7):
    """
    Returns L = elu(y_pred) (considers y_true as 0)
    """

    def true_loss(y_true, y_pred):
        y = -y_pred
        loss = K.elu(y, alpha=alpha)
        return K.sum(loss)

    return true_loss
