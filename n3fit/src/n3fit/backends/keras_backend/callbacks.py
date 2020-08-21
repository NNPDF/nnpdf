"""
    Callbacks to be used during training

    The callbacks defined in this module can be passed to the ``on_epoch_end`` argument
    of the ``perform_fit`` method as a list.
    They must take as input an epoch number and a log of the partial losses
"""

from tensorflow.keras.callbacks import LambdaCallback


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


def gen_stopping_positivity(
    training_model, positivity_datasets, positivity_multipliers, update_freq=100
):
    """
    Updates the given positivity datasets
    with its respective multipliers each ``update_freq`` epochs

    Parameters
    ----------
        training_model: backend Model
            Model being trained
        positivity_datasets: list(str)
            List of the names of the datasets to be trained
        positivity_multipliers: list(float)
            List of multipliers to be applied
        update_freq: int
            each how many epochs the positivity lambda is updated
    """

    if len(positivity_multipliers) != len(positivity_datasets):
        raise ValueError("The number of positivity datasets and multipliers do not match")

    def callback_positivity(epoch, logs):
        if (epoch + 1) % update_freq == 0:
            training_model.multiply_weights(positivity_datasets, positivity_multipliers)

    return LambdaCallback(on_epoch_end=callback_positivity)
