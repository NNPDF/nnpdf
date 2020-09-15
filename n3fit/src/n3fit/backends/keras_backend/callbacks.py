"""
    Callbacks to be used during training

    The callbacks defined in this module can be passed to the ``on_epoch_end`` argument
    of the ``perform_fit`` method as a list.
    They must take as input an epoch number and a log of the partial losses
"""

from tensorflow.keras.callbacks import LambdaCallback, TensorBoard


def gen_stopping_callback(training_model, stopping_object, log_freq=100):
    """
    Given a ``stopping_object``, the callback will monitor the validation chi2
    and will stop the ``training_model`` when the conditions given by ``stopping_object``
    are met.

    Parameters
    ----------
        training_model: backend Model
            Model being trained
        stopping_object: Stopping
            instance of Stopping which controls when the fit should stop
        log_freq: int
            each how manwy epochs the ``print_stats`` argument of ``stopping_object``
            will be set to true
    """
    # TODO: reduce the importance of the callback function moving its logic to Stopping

    def callback_stopping(epoch, logs):
        print_stats = False
        if (epoch + 1) % log_freq == 0:
            print_stats = True
        stopping_object.monitor_chi2(logs, epoch, print_stats=print_stats)
        if stopping_object.stop_here():
            training_model.stop_training = True

    return LambdaCallback(on_epoch_end=callback_stopping)


def gen_lagrange_callback(training_model, datasets, multipliers, update_freq=100):
    """
    Updates the given datasets
    with its respective multipliers each ``update_freq`` epochs

    Parameters
    ----------
        training_model: backend Model
            Model being trained
        datasets: list(str)
            List of the names of the datasets to be trained
        multipliers: list(float)
            List of multipliers to be applied
        update_freq: int
            each how many epochs the positivity lambda is updated
    """

    if len(multipliers) != len(datasets):
        raise ValueError("The number ofvdatasets and multipliers do not match")

    def callback_lagrange(epoch, logs):
        if (epoch + 1) % update_freq == 0:
            training_model.multiply_weights(datasets, multipliers)

    return LambdaCallback(on_epoch_end=callback_lagrange)


def gen_tensorboard_callback(log_dir, profiling=False, histogram_freq = 0):
    """
    Generate tensorboard logging details at ``log_dir``.
    Metrics of the system are saved each epoch.
    If the profiling flag is set to True, it will also attempt
    to save profiling data.

    Note the usage of this callback can hurt performance.

    Parameters
    ----------
        log_dir: str
            Directory in which to save tensorboard details
        profiling: bool
            Whether or not to save profiling information (default False)
    """
    profile_batch = 1 if profiling else 0
    clb = TensorBoard(
        log_dir=log_dir,
        histogram_freq=histogram_freq,
        write_graph=True,
        write_images=False,
        update_freq="epoch",
        profile_batch=profile_batch,
    )
    return clb
