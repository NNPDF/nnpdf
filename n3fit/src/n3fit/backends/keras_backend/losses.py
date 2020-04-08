"""
    Module containing a list of loss functions availables to the fitting code
"""

import tensorflow as tf
from tensorflow.keras import backend as K
from validphys.loader import FallbackLoader as Loader
import numpy as np

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
        if losstype == 'training' and isinstance(spec_dict, dict):
            xgrid_experiment = []
            for dataset_dict in spec_dict.get('datasets'):
                exp_name=dataset_dict.get('name')
                datasetspec = l.check_dataset(exp_name, theoryid=53, cuts="internal") 
                dataset = datasetspec.load()
                xgrid_dataset = dataset.get_kintable()[:, 0]
                xgrid_experiment = np.concatenate((xgrid_experiment, xgrid_dataset))
            trmask = spec_dict.get('trmask')
            xgrid_training =  []
            for num, i in enumerate(trmask.reshape(-1)):
                if i:
                    xgrid_training.append(xgrid_experiment[num])
            xgrid_training = np.array(xgrid_training)
            weight_function = np.sqrt((xgrid_training/5e-6)**0.4)
            tmp = y_true - y_pred
            right_dot = tf.tensordot(invcovmat, K.transpose(tmp/weight_function), axes=1)
            chi2 = tf.tensordot(tmp/weight_function, right_dot, axes=1)
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
