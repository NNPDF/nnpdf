"""
    Library of functions that modify the internal state of Keras/Tensorflow
"""

import os

import psutil

# Despite the current default being tf-eigen, the option below seems to have a positive impact
os.environ.setdefault("KMP_BLOCKTIME", "0")

# Reduce tensorflow verbosity
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "1")
import logging
import random as rn

import keras
from keras import backend as K
import numpy as np

log = logging.getLogger(__name__)

# Prepare Keras-backend dependent functions
if K.backend() in ("torch", "jax"):
    import torch

    def set_eager(flag=True):
        """Pytorch is eager by default"""
        pass

    def set_threading(threads, core):
        """Not implemented"""
        log.info("Setting max number of threads to: %d", threads)
        torch.set_num_threads(threads)

elif K.backend() == "tensorflow":
    import tensorflow as tf

    def set_eager(flag=True):
        """Set eager mode on or off
        for a very slow but fine grained debugging call this function as early as possible
        ideally after the first tf import
        """
        tf.config.run_functions_eagerly(flag)

    def set_threading(threads, cores):
        """Set the Tensorflow inter and intra parallelism options"""
        log.info("Setting the number of cores to: %d", cores)
        try:
            tf.config.threading.set_inter_op_parallelism_threads(threads)
            tf.config.threading.set_intra_op_parallelism_threads(cores)
        except RuntimeError:
            # If tensorflow has already been initiated, the previous calls might fail.
            # This may happen for instance if pdfflow is being used
            log.warning(
                "Could not set tensorflow parallelism settings from n3fit, maybe tensorflow is already initialized by a third program"
            )

else:
    # Keras should've failed by now, if it doesn't it could be a new backend that works ootb?
    log.warning(f"Backend {K.backend()} not recognized. You are entering uncharted territory")


def set_number_of_cores(max_cores=None, max_threads=None):
    """
    Set the maximum number of cores and threads per core to be used by TF.
    It defaults to the number of physical cores
    (and will never surpass it even if max_cores is above)

    Parameters
    ----------
        max_cores: int
            Maximum number of cores to be used
    """
    # Find how many cores we have and how many threads per core
    cores = psutil.cpu_count(logical=False)
    logical = psutil.cpu_count(logical=True)
    tpc = int(logical / cores)

    # We might not have access to all cpus, but assume we get all associated threads for a cpu
    try:
        affinity = psutil.Process().cpu_affinity()
        if len(affinity) != logical:
            cores = int(len(affinity) / tpc)
    except AttributeError:
        # travis Mac OS does not have "cpu_affinity", not sure whether is common to all Macs
        pass

    # In any case, we never want to get above the number provided by the user
    if max_cores is not None:
        cores = min(cores, max_cores)

    threads = tpc * 2
    if max_threads is not None:
        threads = min(max_threads, threads)

    set_threading(threads, cores)


def clear_backend_state():
    """
    Clears the state of the backend.
    Internally it cleans the Keras/TF internal state, liberating the layer names
    and unused memory.
    """
    log.info("Clearing session")
    K.clear_session()


def set_initial_state(debug=False, external_seed=None, max_cores=None, double_precision=False):
    """
    This function sets the initial internal state for the different components of n3fit.

    In debug mode it seeds all seedable libraries, which include:
        - numpy
        - hyperopt
        - python random
        - tensorflow
    The tensorflow/keras part is based on Keras' own
    [guide](https://keras.io/getting_started/faq/#how-can-i-obtain-reproducible-results-using-keras-during-development)
    Note that you might also need PYTHONHASHSEED=0 (outside the program) for full reproducibility.

    To ensure reproducibility in debug mode, if the number of cores is not given,
    it will be set to 1 (with 1 thread per core)

    Parameters
    ----------
        debug: bool
            If this is a debug run, the initial seeds are fixed
        external_seed: int
            Force a seed into numpy, random and tf
        max_cores: int
            Maximum number of cores (as many as physical cores by default)
        double_precision: bool
            If set, use float64 as the default float type
    """
    # If debug mode (or if the external_seed is fixed), fix every non TF seed
    if debug or external_seed is not None:
        if external_seed is None:
            seed = 13
        else:
            seed = external_seed
        log.info("Setting debug seed to: %d", seed)
        # Set the initial seed for the hyperoptimization
        os.environ.setdefault("HYPEROPT_FMIN_SEED", str(seed))

        np.random.seed(seed)
        use_seed = np.random.randint(0, pow(2, 31))
        rn.seed(use_seed)

    # Clear the state of keras in case anyone used it before
    clear_backend_state()

    if double_precision:
        K.set_floatx('float64')

    # Set the number of cores depending on the user choice of max_cores
    # if debug mode and no number of cores set by the user, set to 1
    threads = None  # auto
    if debug and max_cores is None:
        keras.utils.set_random_seed(7331)
        threads = 1
        tf.config.experimental.enable_op_determinism()
    set_number_of_cores(max_cores=max_cores, max_threads=threads)

    # Once again, if in debug mode or external_seed set, set also the TF seed
    if debug or external_seed:
        if K.backend() == "tensorflow":
            # if backend is tensorflow, keep the old seeding
            tf.random.set_seed(use_seed)
        else:
            keras.utils.set_random_seed(use_seed)


def get_physical_gpus():
    """
    Retrieve a list of all physical GPU devices available in the system.

    Returns
    -------
        list: A list of TensorFlow physical devices of type 'GPU'.
    """
    return tf.config.list_physical_devices('GPU')
