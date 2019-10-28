"""
    Library of functions that modify the internal state of Keras/Tensorflow
"""

import numpy as np
import random as rn
import tensorflow as tf
from keras import backend as K


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

    session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=1,
                                inter_op_parallelism_threads=1)
    tf.compat.v1.set_random_seed(use_seed)
    sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)
    #K.set_session(sess)

    return 0


def clear_backend_state():
    """
        Clears the state of the backend and opens a new session.

        Note that this function needs to set the TF session, including threads and processes
        i.e., this function must NEVER be called after setting the initial state.
    """
    print("Clearing session")

    K.clear_session()
    # Don't open threads that you are not going to eat
    session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=2,
                            inter_op_parallelism_threads=8)
    sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)
    #K.set_session(sess)
