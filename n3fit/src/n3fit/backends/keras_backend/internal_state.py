"""
    Library of functions that modify the internal state of Keras/Tensorflow
"""

import os
from psutil import cpu_count

import random as rn
import numpy as np
import tensorflow as tf
from tensorflow.keras import backend as K


def set_initial_state(seed=13):
    """
    Sets the initial state of the backend
    This is the only way of getting reproducible results for keras-tensorflow

    This function needs to be called before __any__ tensorflow related stuff is called so
    it will clear the keras session to ensure the initial state is set

    At the moment this is only enabled for debugging as forces the use of only one thread
    """

    np.random.seed(seed)
    use_seed = np.random.randint(0, pow(2, 31))
    rn.seed(use_seed)

    # Clear the state of keras in case anyone used it before
    K.clear_session()
    tf.config.threading.set_inter_op_parallelism_threads(1)
    tf.config.threading.set_intra_op_parallelism_threads(1)
    tf.random.set_seed(use_seed)

    return 0


def clear_backend_state():
    """
        Clears the state of the backend and opens a new session.

        Note that this function needs to set the TF session, including threads and processes
        i.e., this function must NEVER be called after setting the initial state.
    """
    print("Clearing session")

    # Find how many cores we have and how many threads per core
    cores = cpu_count(logical=False)
    logical = cpu_count(logical=True)
    tpc = int(logical/cores)

    os.environ["KMP_BLOCKTIME"]   = "0"
    os.environ["KMP_SETTINGS"]    = "1"
    os.environ["KMP_AFFINITY"]    = "granularity=fine,verbose,compact,1,0"
    os.environ["OMP_NUM_THREADS"] = str(cores)

    K.clear_session()
    tf.config.threading.set_inter_op_parallelism_threads(tpc)
    tf.config.threading.set_intra_op_parallelism_threads(cores)
