"""
    Callbacks to be used during training

    The callbacks defined in this module can be passed to the ``callbacks`` argument
    of the ``perform_fit`` method as a list.

    For the most typical usage: ``on_epoch_end``,
    they must take as input an epoch number and a log of the partial losses.
"""

import logging
from time import time
import numpy as np
import tensorflow as tf
from tensorflow.keras.callbacks import TensorBoard, Callback

log = logging.getLogger(__name__)


class TimerCallback(Callback):
    """Callback to be used during debugging to time the fit"""

    def __init__(self, count_range=100):
        super().__init__()

        self.all_times = []
        self.every_x = []
        self.x_count = count_range
        self.starting_time = None
        self.last_time = 0

    def on_epoch_end(self, epoch, logs=None):
        """ At the end of every epoch it checks the time """
        new_time = time()
        if epoch == 0:
            # The first epoch is only useful for starting
            self.starting_time = new_time
        else:
            cur_dif = new_time - self.last_time
            self.all_times.append(cur_dif)
            if (epoch + 1) % self.x_count == 0:
                ave = np.mean(self.all_times[-100:])
                log.info(f" > Latest 100 average: {ave:.5} s")
                self.every_x.append(ave)
        self.last_time = new_time

    def on_train_end(self, logs=None):
        """ Print the results """
        total_time = time() - self.starting_time
        n_times = len(self.all_times)
        # Skip the first 100 epochs to avoid fluctuations due to compilations of part of the code
        # by epoch 100 all parts of the code have usually been called so it's a good compromise
        mean = np.mean(self.all_times[min(110, n_times-1):])
        std = np.std(self.all_times[min(110, n_times-1):])
        log.info(f"> > Average time per epoch: {mean:.5} +- {std:.5} s")
        log.info(f"> > > Total time: {total_time/60:.5} min")


class StoppingCallback(Callback):
    """
    Given a ``stopping_object``, the callback will monitor the validation chi2
    and will stop the training model when the conditions given by ``stopping_object``
    are met.

    Parameters
    ----------
        stopping_object: Stopping
            instance of Stopping which controls when the fit should stop
        log_freq: int
            each how many epochs the ``print_stats`` argument of ``stopping_object``
            will be set to true
    """

    def __init__(self, stopping_object, log_freq=100):
        super().__init__()
        self.log_freq = log_freq
        self.stopping_object = stopping_object

    def on_epoch_end(self, epoch, logs=None):
        """ Function to be called at the end of every epoch """
        print_stats = ((epoch + 1) % self.log_freq) == 0
        # Note that the input logs correspond to the fit before the weights are updated
        self.stopping_object.monitor_chi2(logs, epoch, print_stats=print_stats)
        if self.stopping_object.stop_here():
            self.model.stop_training = True

    def on_train_end(self, logs=None):
        """The training can be finished by the stopping or by
        Tensorflow when the number of epochs reaches the maximum.
        In this second case the stopping has to be manually set
        """
        self.stopping_object.make_stop()


class LagrangeCallback(Callback):
    """
    Updates the given datasets
    with its respective multipliers each ``update_freq`` epochs

    Parameters
    ----------
        datasets: list(str)
            List of the names of the datasets to be trained
        multipliers: list(float)
            List of multipliers to be applied
        update_freq: int
            each how many epochs the positivity lambda is updated
    """

    def __init__(self, datasets, multipliers, update_freq=100):
        super().__init__()
        if len(multipliers) != len(datasets):
            raise ValueError("The number of datasets and multipliers do not match")
        self.update_freq = update_freq
        self.datasets = datasets
        self.multipliers = multipliers
        self.updateable_weights = []

    def on_train_begin(self, logs=None):
        """ Save an instance of all relevant layers """
        for layer_name in self.datasets:
            layer = self.model.get_layer(layer_name)
            self.updateable_weights.append(layer.weights)

    @tf.function
    def _update_weights(self):
        """Update all the weight with the corresponding multipliers
        Wrapped with tf.function to compensate the for loops as both weights variables
        and multipliers are known upon first call
        """
        for ws, multiplier in zip(self.updateable_weights, self.multipliers):
            for w in ws:
                w.assign(w * multiplier)

    def on_epoch_end(self, epoch, logs=None):
        """ Function to be called at the end of every epoch """
        if (epoch + 1) % self.update_freq == 0:
            self._update_weights()


def gen_tensorboard_callback(log_dir, profiling=False, histogram_freq=0):
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

class CustomLearningRate(Callback):

    def __init__(self, alphas_settings_dict, epochs):
        power = alphas_settings_dict["learning_rate_settings"]["power"]
        self.decay_steps = alphas_settings_dict["learning_rate_settings"]["decay_steps"]
        initial_learning_rate = alphas_settings_dict["learning_rate_settings"]["initial_learning_rate"]
        end_learning_rate = alphas_settings_dict["learning_rate_settings"]["end_learning_rate"]
        start_fitting_alphas_epoch = alphas_settings_dict["learning_rate_settings"]["start_fitting_alphas_epoch"]
        def decayed_learning_rate(step, decay_steps):
            step = min(step, decay_steps)
            #decay_steps = self.decay_steps * np.ceil(step / self.decay_steps)
            if step < start_fitting_alphas_epoch:
                return 0.0
            else:
                return ((initial_learning_rate - end_learning_rate) *
                        (1 - step / decay_steps) ** (power)
                        ) + end_learning_rate
        self.lr = decayed_learning_rate
        super().__init__()

    def on_epoch_begin(self, epoch, logs=None):
        tf.keras.backend.set_value(
            self.model.optimizer.optimizer_specs[1]["optimizer"].learning_rate,
            self.lr(epoch, self.decay_steps)
            )
