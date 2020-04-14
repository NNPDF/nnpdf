"""
    Module containing a list of loss functions availables to the fitting code
"""

import numpy as np
import tensorflow as tf
from tensorflow.keras import backend as K

from validphys.loader import FallbackLoader as Loader

l = Loader()


DIS_datasets = [
    "NMCPD",
    "NMC",
    "SLACP",
    "SLACD",
    "BCDMSP",
    "BCDMSD",
    "CHORUSNU",
    "CHORUSNB",
    "NTVNUDMN",
    "NTVNBDMN",
    "EMCF2C",
    "HERACOMBNCEM",
    "HERACOMBNCEP460",
    "HERACOMBNCEP575",
    "HERACOMBNCEP820",
    "HERACOMBNCEP920",
    "HERACOMBCCEM",
    "HERACOMBCCEP",
    "HERAF2CHARM",
    "H1HERAF2B",
    "ZEUSHERAF2B",
    "EXTRAPOLATION"
]

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
                if exp_name in DIS_datasets:
                    datasetspec = l.check_dataset(exp_name, theoryid=53, cuts="internal")
                    dataset = datasetspec.load()
                    xgrid_dataset = dataset.get_kintable()[:, 0]
                    xgrid_experiment = np.concatenate((xgrid_experiment, xgrid_dataset))
            if exp_name in DIS_datasets:
                trmask = spec_dict.get("trmask")
                xgrid_training = []
                for num, i in enumerate(trmask.reshape(-1)):
                    if i:
                        xgrid_training.append(xgrid_experiment[num])
                xgrid_training = np.array(xgrid_training)

                def weights(x, polynomial_factors):
                    return sum(polynomial_factors[i] * x ** i for i in range(len(polynomial_factors)))

                def log10(x):
                    numerator = K.log(x)
                    denominator = K.log(tf.constant(10, dtype=numerator.dtype))
                    return numerator / denominator

                DISonly_poly19 = [
                    2.10588812,
                    -4.88688822,
                    -40.4326353,
                    -220.741839,
                    -681.465768,
                    -1297.64601,
                    -1654.86188,
                    -1488.34315,
                    -972.269594,
                    -467.095599,
                    -164.20549,
                    -40.8970087,
                    -6.49081327,
                    -0.35816766,
                    0.110690877,
                    0.0345086403,
                    0.00494895461,
                    0.000416682047,
                    1.98498358e-05,
                    4.1579982e-07,
                ]

                weight_function = 10 ** 0.5 * (10 ** weights(log10(xgrid_training), DISonly_poly19)) ** (
                    -0.2
                )

                tmp = y_true - y_pred
                weight_function = tf.cast(weight_function, dtype=tmp.dtype)
                tmp *= weight_function

                right_dot = tf.tensordot(invcovmat, K.transpose(tmp), axes=1)
                chi2 = tf.tensordot(tmp, right_dot, axes=1)
            else:
                tmp = y_true - y_pred
                right_dot = tf.tensordot(invcovmat, K.transpose(tmp), axes=1)
                chi2 = tf.tensordot(tmp, right_dot, axes=1)
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
