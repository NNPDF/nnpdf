"""
    Module containing a list of loss functions availables to the fitting code
"""

import numpy as np
import tensorflow as tf
from tensorflow.keras import backend as K

from validphys.loader import FallbackLoader as Loader

l = Loader()


def l_invcovmat(invcovmat_np, losstype="validation", exp_name=None, spec_dict=None):
    """
    Returns a loss function such that:
    L = \sum_{ij} (yt - yp)_{i} invcovmat_{ij} (yt - yp)_{j}
    """
    invcovmat = K.constant(invcovmat_np)

    def true_loss(y_true, y_pred):
        # (yt - yp) * covmat * (yt - yp)
        # NOTE: if isinstance isn't here we get a NameError: trainmask is not defined
        if losstype == "training" and isinstance(spec_dict, dict):
            xgrid_experiment = []
            for dataset_dict in spec_dict.get("datasets"):
                exp_name = dataset_dict.get("name")
                datasetspec = l.check_dataset(exp_name, theoryid=53, cuts="internal")
                dataset = datasetspec.load()
                xgrid_dataset = dataset.get_kintable()[:, 0]
                xgrid_experiment = np.concatenate((xgrid_experiment, xgrid_dataset))
            trmask = spec_dict.get("trmask")
            xgrid_training = []
            for num, i in enumerate(trmask.reshape(-1)):
                if i:
                    xgrid_training.append(xgrid_experiment[num])
            xgrid_training = np.array(xgrid_training)
            # weight_function = np.sqrt((xgrid_training/5e-6)**0.4)
            def weights(x):
                polynomial_density_fit = (
                    2.1057592888948253
                    + -4.164527162015765 * x
                    + -26.400855363117646 * x ** 2
                    + -122.75480240246661 * x ** 3
                    + -327.1547368943516 * x ** 4
                    + -523.4899128971089 * x ** 5
                    + -543.0024185862545 * x ** 6
                    + -384.4247519539939 * x ** 7
                    + -191.61198802145452 * x ** 8
                    + -68.34473527884417 * x ** 9
                    + -17.514055465376153 * x ** 10
                    + -3.1967717745260646 * x ** 11
                    + -0.4053698357109279 * x ** 12
                    + -0.033923870763836385 * x ** 13
                    + -0.00168403639998969 * x ** 14
                    + -3.754716120287798e-05 * x ** 15
                )
                return polynomial_density_fit
            def log10(x):
                numerator = K.log(x)
                denominator = K.log( tf.constant(10, dtype=numerator.dtype))
                return numerator/denominator
            weight_function = (10**weights(log10(xgrid_training)))**(-0.5)
            
            tmp = y_true - y_pred
            weight_function = tf.cast(weight_function,dtype=tmp.dtype)

            right_dot = tf.tensordot(invcovmat, K.transpose(tmp * weight_function), axes=1)
            chi2 = tf.tensordot(tmp * weight_function, right_dot, axes=1)
        else:
            tmp = y_true - y_pred
            right_dot = tf.tensordot(invcovmat, K.transpose(tmp), axes=1)
            chi2 = tf.tensordot(tmp, right_dot, axes=1)
        return chi2

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

def l_diaginvcovmat(diaginvcovmat_np):
    """
    Returns a loss function such that:
    L = sum_{i} (yt - yp)_{i} invcovmat_{ii} (yt - yp)_{i}
    diaginvcovmat_np should be 1d
    """
    invcovmat = K.constant(diaginvcovmat_np)

    def true_loss(y_true, y_pred):
        tmp = y_true - y_pred
        return tf.tensordot(invcovmat, K.transpose(tmp*tmp), axes=1)
    return true_loss
