def set_initial_state(debug = False, seed = 13):
    """
    Setting the initial state for debugging can be very tricky
    This is the only way of getting reproducible results for keras-tensorflow
    Note that this function needs to be called before _any_ tensorflow related stuff is called

    At the moment this is only enabled for debugging as forces the use of only one thread
    """
    if not debug:
        return 0

    import numpy as np
    import random as rn

    np.random.seed(seed)
    use_seed = np.random.randint(0, pow(2,31))
    rn.seed(use_seed)

    from keras import backend as K
    session_conf = K.tf.ConfigProto(intra_op_parallelism_threads=1,
                                inter_op_parallelism_threads=1)
    K.tf.set_random_seed(use_seed)
    sess = K.tf.Session(graph=K.tf.get_default_graph(), config=session_conf)
    K.set_session(sess)

    return 0

def clear_backend_state(debug):
    if not debug:
        print("Clearing session")
        from keras import backend as K
        K.clear_session()
        # Don't open threads that you are not going to eat
        session_conf = K.tf.ConfigProto(intra_op_parallelism_threads=2,
                                inter_op_parallelism_threads=8)
        sess = K.tf.Session(graph=K.tf.get_default_graph(), config=session_conf)
        K.set_session(sess)
        
