"""
Callbacks to be used during training

The callbacks defined in this module can be passed to the ``callbacks`` argument
of the ``perform_fit`` method as a list.

For the most typical usage: ``on_batch_end``,
they must take as input an epoch number and a log of the partial losses.

Note: the terminology used everywhere refers to a single training step as a single epoch.
It turns out that to avoid tensorflow overhead, it is beneficial to write a step as a
single batch instead. So callbacks must use ``on_batch_end``.
"""

import logging
from time import time

from keras import backend as K
from keras.callbacks import Callback, TensorBoard
import numpy as np
import tensorflow as tf
import pandas as pd

from .operations import decorator_compiler

log = logging.getLogger(__name__)


class CallbackStep(Callback):
    """
    Wrapper around the keras Callback that keeps track of how the steps are divided
    between epochs and batches.
    The callback will call ``on_step_end`` instead of ``on_batch_end``.
    """

    def __init__(self):
        super().__init__()
        self.steps_in_epoch = 0
        self.epochs_finished = 0
        self.steps_per_epoch = 0  # will be defined in the first epoch
        self._previous_logs = {}

    def on_epoch_end(self, epoch, logs=None):
        if self.steps_per_epoch == 0:
            self.steps_per_epoch = self.steps_in_epoch
        self.steps_in_epoch = 0
        self.epochs_finished += 1

    def on_batch_end(self, batch, logs=None):
        #tf.print("CallbackStep.on_batch_end called, batch=", batch)
        #if logs is not None:
        #    tf.print("logs keys:", list(logs.keys()))
        step_number = self.steps_in_epoch + self.epochs_finished * self.steps_per_epoch
        self.on_step_end(step_number, logs)
        self.steps_in_epoch += 1

    def correct_logs(self, logs: dict) -> dict:
        """
        The logs that get computed by default are an average over batches.
        This converts it into the logs for the current step.
        """
        #tf.print("correct_logs called, logs keys:", list(logs.keys()) if logs else 'None')
        corrected_logs = {}
        for k in logs.keys():
            previous_total = self._previous_logs.get(k, 0.0) * self.steps_in_epoch
            current_total = logs[k] * (self.steps_in_epoch + 1)
            corrected_logs[k] = current_total - previous_total
        self._previous_logs = logs.copy()
        return corrected_logs


class TimerCallback(CallbackStep):
    """Callback to be used during debugging to time the fit"""

    def __init__(self, count_range=100):
        super().__init__()

        self.all_times = []
        self.every_x = []
        self.x_count = count_range
        self.starting_time = None
        self.last_time = 0

    def on_step_end(self, epoch, logs=None):
        """At the end of every epoch it checks the time"""
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
        """Print the results"""
        total_time = time() - self.starting_time
        n_times = len(self.all_times)
        # Skip the first 100 epochs to avoid fluctuations due to compilations of part of the code
        # by epoch 100 all parts of the code have usually been called so it's a good compromise
        mean = np.mean(self.all_times[min(110, n_times - 1) :])
        std = np.std(self.all_times[min(110, n_times - 1) :])
        log.info(f"> > Average time per epoch: {mean:.5} +- {std:.5} s")
        log.info(f"> > > Total time: {total_time/60:.5} min")


class StoppingCallback(CallbackStep):
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

    def on_step_end(self, epoch, logs=None):
        """Function to be called at the end of every epoch
        Every ``log_freq`` number of epochs, the ``monitor_chi2`` method of the ``stopping_object``
        will be called and the validation loss (broken down by experiment) will be logged.
        For the training model only the total loss is logged during the training.
        """
        print_stats = ((epoch + 1) % self.log_freq) == 0
        # Note that the input logs correspond to the fit before the weights are updated
        logs = self.correct_logs(logs)

        # WARNING: this line seems to be necessary for jax
        # otherwise the validation model itself cannot run compute_losses
        # but it needs to be run every epoch, which makes no sense
        if K.backend() == "jax":
            _ = self.model.compute_losses()

        self.stopping_object.monitor_chi2(logs, epoch, print_stats=print_stats)
        if self.stopping_object.stop_here():
            self.model.stop_training = True

    def on_train_end(self, logs=None):
        """The training can be finished by the stopping or by
        Tensorflow when the number of epochs reaches the maximum.
        In this second case the stopping has to be manually set
        """
        self.stopping_object.make_stop()


class LagrangeCallback(CallbackStep):
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
        """Save an instance of all relevant layers"""
        for layer_name in self.datasets:
            layer = self.model.get_layer(layer_name)
            self.updateable_weights.append(layer.weights)

    @decorator_compiler
    def _update_weights(self):
        """Update all the weight with the corresponding multipliers
        Wrapped with tf.function to compensate the for loops as both weights variables
        and multipliers are known upon first call
        """
        for ws, multiplier in zip(self.updateable_weights, self.multipliers):
            for w in ws:
                w.assign(w * multiplier)

    def on_step_end(self, epoch, logs=None):
        """Function to be called at the end of every epoch"""
        if (epoch + 1) % self.update_freq == 0:
            self._update_weights()

class KLAnnealingCallback(CallbackStep):
    """
    Updates the tensorflow variable "kl_beta" to increase 
    from 0 to 1 over a specified number of warmup steps.

    Parameters
    ----------
    kl_beta: tf.Variable
        The variable used in the loss function to scale the KL term.
    warmup_steps: int
        The number of steps (epochs/batches) until kl_beta reaches 1.0.
    """

    def __init__(self, kl_beta, warmup_steps=2000, replica_path=None):
        super().__init__()
        self.kl_beta = kl_beta
        self.warmup_steps = max(1, warmup_steps)
        self.total_steps = tf.Variable(0, trainable=False, dtype=tf.float32)
        self.klarray = []
        self.replica_path = replica_path


    def linear_annealing(self, step, logs=None):
        """
        Updates the beta variable at the end of each epoch
        Tracks total steps directly to bypass correct_logs issues.
        """
        # Linearly increase beta until it hits 1.0
        self.new_beta = min(1.0, self.total_steps / self.warmup_steps)
        self.kl_beta.assign(tf.cast(self.new_beta, self.kl_beta.dtype))
        if self.total_steps % 100 == 0:  # debug log every 500 steps
            tf.print("[INFO]: Linear KLAnnealing: step=", self.total_steps, " beta=", self.new_beta)

    def sigmoid_annealing(self, step, logs=None):
        k = 0.0005
        t0 = self.warmup_steps / 2
        
        self.new_beta = 1 / (1 + tf.exp(-k*(self.total_steps-t0)))
        self.kl_beta.assign(tf.cast(self.new_beta, self.kl_beta.dtype))
        if self.total_steps % 100 == 0:  # debug log every 500 steps
            tf.print("[INFO]: Sigmoid KLAnnealing: step=", self.total_steps, " beta=", self.new_beta)

    def cyclical_annealing(self, step, logs=None):
        
        cycle_len = self.warmup_steps // 4
        # Linear ramp within the cycle
        rel_step = self.total_steps % cycle_len
        self.new_beta = min(1.0, rel_step / (cycle_len * 0.5)) # Ramp for 50% of the cycle
        self.kl_beta.assign(tf.cast(self.new_beta, self.kl_beta.dtype))
        if self.total_steps % 100 == 0:  # debug log every 500 steps
            tf.print("[INFO]: Cyclical KLAnnealing: step=", self.total_steps, " beta=", self.new_beta)
    
    def old_cosine_annealing(self, step, logs=None):
        beta_min = tf.cast(0.1, tf.float32)
        beta_max = tf.cast(1.0, tf.float32)
        ratio = tf.cast(0.5, tf.float32)
        cycle_len = tf.cast(self.warmup_steps // 4, tf.float32)
        
        t_mod = tf.cast(self.total_steps % cycle_len, tf.float32)
        ramp_end = ratio * cycle_len
        
        phi = t_mod / ramp_end
        ramped = beta_min + 0.5 * (beta_max - beta_min) * (1 - tf.math.cos(np.pi * phi))
        
        self.new_beta = tf.cond(t_mod < ramp_end, lambda: ramped, lambda: beta_max)   

    def cosine_annealing(self, step, logs=None):
        beta_min = tf.cast(0.1, tf.float32)
        beta_max = tf.cast(1.0, tf.float32)
        ratio = tf.cast(0.5, tf.float32)
        cycle_len = tf.cast(self.warmup_steps // 4, tf.float32)
        
        t_mod = tf.cast(self.total_steps % cycle_len, tf.float32)
        ramp_end = ratio * cycle_len
        
        phi = t_mod / ramp_end
        ramped = beta_min + 0.5 * (beta_max - beta_min) * (1 - tf.math.cos(np.pi * phi))
        
        cyclical_beta = tf.cond(t_mod < ramp_end, lambda: ramped, lambda: beta_max)

        self.new_beta = tf.cond(step > self.warmup_steps, lambda: beta_max, lambda: cyclical_beta)
        
    def on_step_end(self, step, logs=None):
        self.total_steps.assign_add(1)
        #tf.print("KL.on_step_end: total_steps=", self.total_steps, " model=", self.model is not None)
        self.cosine_annealing(step)
        old_beta = self.kl_beta.read_value()
        self.kl_beta.assign(tf.cast(self.new_beta, self.kl_beta.dtype))
        new_beta = self.kl_beta.read_value()
        self.klarray.append(new_beta)
    
    def on_train_end(self, logs=None):
        path = self.replica_path / "kl.csv" if self.replica_path is not None else "kl.csv"
        pd.DataFrame(self.klarray, columns=["kl_beta"]).to_csv(path, index=False)


def gen_tensorboard_callback(log_dir, profiling=False, histogram_freq=0):
    """
    Generate tensorboard logging details at ``log_dir``.
    Metrics of the system are saved each epoch.
    If the profiling flag is set to True, it will also attempt
    to save profiling data.

    Note the usage of this callback can hurt performance
    At the moment can only be used with TensorFlow: https://github.com/keras-team/keras/issues/19121

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
