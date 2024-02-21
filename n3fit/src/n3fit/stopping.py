"""
    Module containing the classes related to the stopping alogirthm

    In this module there are four Classes:

    - FitState: this class contains the information of the fit
            for a given point in history
    - FitHistory: this class contains the information necessary
            in order to reset the state of the fit to the point
            in which the history was saved.
            i.e., a list of FitStates
    - Stopping: this class monitors the chi2 of the validation
            and training sets and decides when to stop
    - Positivity: Decides whether a given point fullfills the positivity conditions
    - Validation: Controls the NNPDF cross-validation algorithm

    Note:
        There are situations in which the validation set is empty, in those cases
    the training set is used as validation set.
    This implies several changes in the behaviour of this class as the training chi2 will
    now be monitored for stability.
        In order to parse the set of loss functions coming from the backend::MetaModel,
    the function `normalize_chi2s` relies on the fact that they are all suffixed with `_loss`
    the validation case, instead, is suffixed with `val_loss`. In the particular casse in
    which both training and validation model correspond to the same backend::MetaModel only
    the `_loss` suffix can be found. This is taken into account by the class `Stopping`
    which will tell `Validation` that no validation set was found and that the training is to
    be used instead.
"""
import logging

import numpy as np

from n3fit.backends import LossMetric

log = logging.getLogger(__name__)

# Put a very big number here so that we for sure discard this run
# AND we have a clear marker that something went wrong, not just a bad fit
TERRIBLE_CHI2 = 1e10
INITIAL_CHI2 = 1e9

# Pass/veto keys
POS_OK = "POS_PASS"
POS_BAD = "POS_VETO"
THRESHOLD_POS = 1e-6


def normalize_chi2s(loss_dict, chi2_point_dict, suffix="loss"):
    """
    From a dict of all losses, select those that contribute to the chi2 (i.e. exclude positivity
    and integrability penalties).
    Normalize by the number of points in each dataset.

    Parameters
    ----------
        loss_dict: dict
            dictionary of losses per experiment
        chi2_point_dict: dict
            dictionary with the name of the experiments to be taken into account
            and the number of datapoints of the experiments
        suffix: str (default: ``loss``)
            suffix of the loss layer, Keras default is _loss

    Returns
    -------
        dict_chi2: dict
            dictionary of {'expname': loss }, including the normalized total chi2 as 'total',
                and the loss including penalties as 'loss'
    """
    dict_chi2 = {}
    total_points = 0
    total_loss = 0
    for exp_name, npoints in chi2_point_dict.items():
        loss = np.array(loss_dict[exp_name + f"_{suffix}"])
        dict_chi2[exp_name] = loss / npoints
        total_points += npoints
        total_loss += loss

    dict_chi2["total"] = total_loss / total_points
    dict_chi2["loss"] = loss_dict["loss"]

    return dict_chi2


class FitState:
    """
    Holds the state of the chi2 during the fit, for all replicas and one epoch

    Note: the training chi2 is computed before the update of the weights
    so it is the chi2 that informed the updated corresponding to this state.
    The validation chi2 instead is computed after the update of the weights.

    Parameters
    ----------
        training_chi2s: dict
            chi2 losses for the training model
        validation_chi2s: dict
            chi2 losses for the validation model
    """

    def __init__(self, training_chi2s, validation_chi2s):
        # these are summed over replicas
        self.all_tr_chi2 = training_chi2s
        # these are per replica
        self.all_vl_chi2 = validation_chi2s

        self.tr_chi2 = self.all_tr_chi2["total"]
        self.vl_chi2 = self.all_vl_chi2["total"]

        self.tr_loss = self.all_tr_chi2["loss"]
        self.vl_loss = self.all_vl_chi2["loss"]
        del self.all_tr_chi2["loss"]
        del self.all_vl_chi2["loss"]

    def all_tr_chi2_for_replica(self, i_replica):
        """Return the tr chi2 per dataset for a given replica"""
        return {k: np.take(v, i_replica) for k, v in self.all_tr_chi2.items()}

    def all_vl_chi2_for_replica(self, i_replica):
        """Return the vl chi2 per dataset for a given replica"""
        return {k: np.take(v, i_replica) for k, v in self.all_vl_chi2.items()}

    def total_partial_tr_chi2(self):
        """Return the tr chi2 summed over replicas per experiment"""
        return {k: np.sum(v) for k, v in self.all_tr_chi2.items()}

    def total_partial_vl_chi2(self):
        """Return the vl chi2 summed over replicas per experiment"""
        return {k: np.sum(v) for k, v in self.all_tr_chi2.items()}

    def total_tr_chi2(self):
        """Return the total tr chi2 summed over replicas"""
        return np.sum(self.tr_chi2)

    def total_vl_chi2(self):
        """Return the total vl chi2 summed over replicas"""
        return np.sum(self.vl_chi2)

    def __str__(self):
        return f"chi2: tr={self.tr_chi2} vl={self.vl_chi2}"


class FitHistory:
    """
    Keeps a list of FitState items holding the full chi2 history of the fit.

    Parameters
    ----------
        tr_ndata: dict
            dictionary of {dataset: n_points} for the training data
        vl_ndata: dict
            dictionary of {dataset: n_points} for the validation data
    """

    def __init__(self):
        # Save a list of status for the entire fit
        self._history = []
        self.final_epoch = None

    def get_state(self, epoch):
        """Get the FitState of the system for a given epoch"""
        try:
            return self._history[epoch]
        except IndexError as e:
            raise ValueError(
                f"Tried to get obtain the state for epoch {epoch} when only {len(self._history)} epochs have been saved"
            ) from e

    def register(self, epoch, training_chi2s, validation_chi2s):
        """Save a new fitstate and updates the current final epoch

        Parameters
        ----------
            epoch: int
                the current epoch of the fit
            training_chi2s: dict
                chi2 losses for the training model
            validation_chi2s: dict
                chi2 losses for the validation model

        Returns
        -------
            FitState
        """
        # Save all the information in a fitstate object
        fitstate = FitState(training_chi2s, validation_chi2s)
        self.final_epoch = epoch
        self._history.append(fitstate)
        return fitstate


class Stopping:
    """
    Driver of the stopping algorithm

    Note, if the total number of points in the validation dictionary is None, it is assumed
    the validation_model actually corresponds to the training model.

    Parameters
    ----------
        observables_model: n3fit.backends.MetaModel
            The model outputting the observables
        training_losses: dict
            dictionary of {dataset: loss} for the training data
        all_data_dicts: dict
            tuple of dicts containing info about training, validation, stopping datasets
        pdf_model: n3fit.backends.MetaModel
           pdf_model being trained
        threshold_positivity: float
           maximum value allowed for the sum of all positivity losses
        total_epochs: int
           total number of epochs
        stopping_patience: int
           how many epochs to wait for the validation loss to improve
        threshold_chi2: float
            maximum value allowed for chi2
        dont_stop: bool
           dont care about early stopping
    """

    def __init__(
        self,
        observables_model,
        training_losses,
        all_data_dicts,
        pdf_model,
        threshold_positivity=THRESHOLD_POS,
        total_epochs=0,
        stopping_patience=7000,
        threshold_chi2=10.0,
        dont_stop=False,
    ):
        self._observables_model = observables_model
        self._pdf_model = pdf_model

        # Save the validation object
        self._training_losses = training_losses

        # Create the History object
        tr_ndata, vl_ndata, pos_sets = all_data_dicts
        self.tr_ndata = tr_ndata
        self.vl_ndata = vl_ndata if vl_ndata is not None else tr_ndata
        self.tr_suffix = "loss"
        self.vl_suffix = "val_loss" if vl_ndata is not None else "loss"
        self._history = FitHistory()

        # And the positivity checker
        self._positivity = Positivity(threshold_positivity, pos_sets)

        # Initialize internal variables for the stopping
        self._n_replicas = pdf_model.num_replicas
        self._threshold_chi2 = threshold_chi2
        self._stopping_degrees = np.zeros(self._n_replicas, dtype=int)
        self._counts = np.zeros(self._n_replicas, dtype=int)

        self._dont_stop = dont_stop
        self._stop_now = False
        self.stopping_patience = stopping_patience
        self.total_epochs = total_epochs

        self._stop_epochs = [total_epochs - 1] * self._n_replicas
        self._best_epochs = [-1] * self._n_replicas
        self.positivity_statuses = [POS_BAD] * self._n_replicas
        self._best_weights = [None] * self._n_replicas
        self._best_val_chi2s = [INITIAL_CHI2] * self._n_replicas

        self.previous_weights = self._pdf_model.get_all_replica_weights()

    @property
    def vl_chi2(self):
        """Current validation chi2"""
        validation_losses = self.compute_validation_losses()
        validation_chi2s = normalize_chi2(validation_losses, self.vl_ndata, self.vl_suffix)
        return validation_chi2s["total"]

    @property
    def e_best_chi2(self):
        """Epoch of the best chi2, if there is no best epoch, return last"""
        best_or_last_epochs = [
            best if best is not None else last
            for best, last in zip(self._best_epochs, self._stop_epochs)
        ]
        return best_or_last_epochs

    @property
    def stop_epoch(self):
        """Epoch in which the fit is stopped"""
        return -1 if self._history.final_epoch is None else self._history.final_epoch + 1

    @property
    def positivity_status(self):
        """Returns POS_PASS if positivity passes or veto if it doesn't
        for each replica"""
        return self.positivity_statuses

    def compute_validation_losses(self):
        """
        Returns dict {loss_name: per replica losses}

        Used inside training loop, no observables_model call done
        """
        metrics = {
            m.name: m.per_replica_losses.numpy()
            for m in self._observables_model.metrics
            if isinstance(m, LossMetric)
        }
        metrics["loss"] = np.sum([m for m in metrics.values()], axis=0)

        # TODO Aron: temporary fix for some issue with POS suffix.
        metrics2 = {}
        for name in metrics.keys():
            if name[:3] == 'POS':
                newname = name.replace('val', 'tr')
            else:
                newname = name
            metrics2[newname] = metrics[name]

        metrics = metrics2
        return metrics

    def monitor_chi2(self, training_losses, epoch, print_stats=False):
        """
        Function to be called at the end of every epoch.
        Stores the total chi2 of the training set as well as the
        total chi2 of the validation set.
        If the training chi2 is below a certain threshold,
        stores the state of the model which gave the minimum chi2
        as well as the epoch in which occurred
        If the epoch is a multiple of save_all_each then we also save the per-exp chi2

        Returns True if the run seems ok and False if a NaN is found

        Parameters
        ----------
            training_losses: dict
                output of a .fit() call, dictionary of the total loss (summed over replicas) for
                each experiment
            epoch: int
                index of the epoch

        Returns
        -------
            pass_ok: bool
                true/false according to the status of the run
        """
        # Step 1. Check whether the fit has NaN'd and stop it if so
        if np.isnan(training_losses["loss"]):
            log.warning(" > NaN found, stopping activated")
            self.make_stop()
            return False

        # Step 2. Compute the validation metrics
        validation_losses = self.compute_validation_losses()

        # Step 3. Register the current point in (the) history
        training_chi2s = normalize_chi2s(training_losses, self.tr_ndata, self.tr_suffix)
        validation_chi2s = normalize_chi2s(validation_losses, self.vl_ndata, self.vl_suffix)
        fitstate = self._history.register(epoch, training_chi2s, validation_chi2s)
        if print_stats:
            self.print_current_stats(epoch, fitstate)

        # Step 4. Check whether this is a better fit
        #         this means improving vl_chi2 and passing positivity
        # Don't start counting until the chi2 of the validation goes below a certain threshold
        # once we start counting, don't bother anymore
        passes = self._counts | (fitstate.vl_chi2 < self._threshold_chi2)
        passes &= fitstate.vl_loss < self._best_val_chi2s
        # And the ones that pass positivity
        passes &= self._positivity(validation_losses)

        self._stopping_degrees += self._counts

        # Step 5. loop over the valid indices to check whether the vl improved
        for i_replica in np.where(passes)[0]:
            self._best_epochs[i_replica] = epoch - 1  # validation is from previous epoch
            # By definition, if we have a ``best_epoch`` then positivity passed
            self.positivity_statuses[i_replica] = POS_OK

            self._best_weights[i_replica] = self.previous_weights[i_replica]
            self._best_val_chi2s[i_replica] = self._history.get_state(epoch).vl_loss[i_replica]

            self._stopping_degrees[i_replica] = 0
            self._counts[i_replica] = 1

        stop_replicas = self._counts & (self._stopping_degrees > self.stopping_patience)
        for i_replica in np.where(stop_replicas)[0]:
            self._stop_epochs[i_replica] = epoch
            self._counts[i_replica] = 0

        # By using the stopping degree we only stop when none of the replicas are improving anymore
        if min(self._stopping_degrees) > self.stopping_patience:
            self.make_stop()

        self.previous_weights = self._pdf_model.get_all_replica_weights()

        return True

    def make_stop(self):
        """Convenience method to set the stop_now flag
        and reload the history to the point of the best model if any
        """
        self._stop_now = True
        self._restore_best_weights()

    def _restore_best_weights(self):
        for i_replica, weights in enumerate(self._best_weights):
            if weights is not None:
                self._pdf_model.set_replica_weights(weights, i_replica)

    def print_current_stats(self, epoch, fitstate):
        """
        Prints ``fitstate`` training and validation chi2s
        """
        epoch_index = epoch + 1
        tr_chi2 = fitstate.total_tr_chi2()
        vl_chi2 = fitstate.total_vl_chi2()
        total_str = f"At epoch {epoch_index}/{self.total_epochs}, total chi2: {tr_chi2}\n"

        # The partial chi2 makes no sense for more than one replica at once:
        if self._n_replicas == 1:
            partial_tr_chi2 = fitstate.total_partial_tr_chi2()
            partials = []
            for experiment, chi2 in partial_tr_chi2.items():
                partials.append(f"{experiment}: {chi2:.3f}")
            total_str += ", ".join(partials) + "\n"
        total_str += f"Validation chi2 at this point: {vl_chi2}"
        log.info(total_str)

    def stop_here(self):
        """Returns the stopping status
        If `dont_stop` is set returns always False (i.e., never stop)
        """
        if self._dont_stop:
            return False
        else:
            return self._stop_now

    def chi2exps_json(self, i_replica=0, log_each=100):
        """
        Returns and apt-for-json dictionary with the status of the fit every `log_each` epochs

        Parameters
        ----------
            i_replica: int
                which replica are we writing the log for
            log_each: int
                every how many epochs to print the log

        Returns
        -------
            file_list: list(str)
                a list of strings to be printed as `chi2exps.log`
        """
        final_epoch = self._history.final_epoch
        json_dict = {}

        for epoch in range(log_each - 1, final_epoch + 1, log_each):
            fitstate = self._history.get_state(epoch)
            all_tr = fitstate.all_tr_chi2_for_replica(i_replica)
            all_vl = fitstate.all_vl_chi2_for_replica(i_replica)

            tmp = {exp: {"training": tr_chi2} for exp, tr_chi2 in all_tr.items()}
            for exp, vl_chi2 in all_vl.items():
                if exp not in tmp:
                    tmp[exp] = {"training": None}
                tmp[exp]["validation"] = vl_chi2

            json_dict[epoch + 1] = tmp
        return json_dict


class Positivity:
    """
    Controls the positivity requirements.

    In order to check the positivity passes will check the history of the fitting
    as the fitting included positivity sets.
    If the sum of all positivity sets losses is above a certain value the model is
    not accepted and the training continues.

    Parameters
    ----------
        threshold_positivity: float
            maximum value allowed for the sum of all positivity losses
        positivity_sets: list
            list of positivity datasets
    """

    def __init__(self, threshold, positivity_sets):
        self.threshold = threshold
        self.positivity_sets = positivity_sets

    def check_positivity(self, history_object):
        """
                This function receives a history objects and loops over the
                positivity_sets to check the value of the positivity loss.

                If the positivity loss is above the threshold, the positivity fails
                otherwise, it passes.
                It returns an array booleans which are True if positivity passed
        story_object[key_loss] < self.threshold
                Parameters
                ----------
                    history_object: dict
                        dictionary of entries in the form  {'name': loss}, output of a MetaModel .fit()
        """
        positivity_pass = True
        for key in self.positivity_sets:
            key_loss = f"{key}_loss"
            positivity_pass &= history_object[key_loss] < self.threshold
        return np.array(positivity_pass)

    def __call__(self, validation_losses):
        """
        Checks whether a given FitState object
        passes the positivity requirement
        """
        return self.check_positivity(validation_losses)
