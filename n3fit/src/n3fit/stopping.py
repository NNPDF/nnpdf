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
    the function `parse_losses` relies on the fact that they are all suffixed with `_loss`
    the validation case, instead, is suffixed with `val_loss`. In the particular casse in
    which both training and validation model correspond to the same backend::MetaModel only
    the `_loss` suffix can be found. This is taken into account by the class `Stopping`
    which will tell `Validation` that no validation set was found and that the training is to
    be used instead.
"""
import logging

import numpy as np

log = logging.getLogger(__name__)

# Put a very big number here so that we for sure discard this run
# AND we have a clear marker that something went wrong, not just a bad fit
TERRIBLE_CHI2 = 1e10
INITIAL_CHI2 = 1e9

# Pass/veto keys
POS_OK = "POS_PASS"
POS_BAD = "POS_VETO"
THRESHOLD_POS = 1e-6


def parse_ndata(all_data):
    """
    Parses the list of dictionaries received from ModelTrainer
    into a dictionary containing only the name of the experiments
    together with the number of points.

    Returns
    -------
        `tr_ndata`
            dictionary of {'exp' : ndata}
        `vl_ndata`
            dictionary of {'exp' : ndata}
        `pos_set`: list of the names of the positivity sets

    Note: if there is no validation (total number of val points == 0)
    then vl_ndata will point to tr_ndata
    """
    tr_ndata_dict = {}
    vl_ndata_dict = {}
    pos_set = []
    for dictionary in all_data:
        exp_name = dictionary["name"]
        if dictionary.get("count_chi2"):
            tr_ndata = dictionary["ndata"]
            vl_ndata = dictionary["ndata_vl"]
            if tr_ndata:
                tr_ndata_dict[exp_name] = tr_ndata
            if vl_ndata:
                vl_ndata_dict[exp_name] = vl_ndata
        if dictionary.get("positivity") and not dictionary.get("integrability"):
            pos_set.append(exp_name)
    if not vl_ndata_dict:
        vl_ndata_dict = None
    return tr_ndata_dict, vl_ndata_dict, pos_set


def parse_losses(history_object, data, suffix="loss"):
    """
    Receives an object containing the chi2
    Usually a history object, but it can come in the form of a dictionary.

    It loops over the dictionary and uses the npoints_data dictionary to
    normalize the chi2 and return backs a tuple (`total`, `tr_chi2`)

    Parameters
    ----------
        history_object: dict
            A history object dictionary
        data: dict
            dictionary with the name of the experiments to be taken into account
            and the number of datapoints of the experiments
        suffix: str (default: ``loss``)
            suffix of the loss layer, Keras default is _loss

    Returns
    -------
        total_loss: float
            Total value for the loss
        dict_chi2: dict
            dictionary of {'expname' : loss }
    """
    try:
        hobj = history_object.history
    except AttributeError:  # So it works whether we pass the out or the out.history
        hobj = history_object

    # In the general case epochs = 1.
    # In case that we are doing more than 1 epoch, take the last result as it is the result
    # the model is in at the moment
    # This value is only used for printing output purposes so should not have any significance
    dict_chi2 = {}
    total_points = 0
    total_loss = 0
    for exp_name, npoints in data.items():
        loss = np.array(hobj[exp_name + f"_{suffix}"])
        dict_chi2[exp_name] = loss / npoints
        total_points += npoints
        total_loss += loss

    # By taking the loss from the history object we would be saving the total loss
    # including positivity sets and (if added/enabled) regularizsers
    # instead we want to restrict ourselves to the loss coming from experiments
    # total_loss = np.mean(hobj["loss"]) / total_points
    total_loss /= total_points
    dict_chi2["total"] = total_loss
    return total_loss, dict_chi2


class FitState:
    """
    Holds the state of the chi2 during the fit for all replicas

    It holds the necessary information to reload the fit
    to a specific point in time if we are interested on reloading
    (otherwise the relevant variables stay empty to save memory)

    Note: the training chi2 is computed before the update of the weights
    so it is the chi2 that informed the updated corresponding to this state.
    The validation chi2 instead is computed after the update of the weights.

    Parameters
    ----------
        training_info: dict
            all losses for the training model
        validation_info: dict
            all losses for the validation model
    """

    vl_ndata = None
    tr_ndata = None
    vl_suffix = None

    def __init__(self, training_info, validation_info):
        if self.vl_ndata is None or self.tr_ndata is None or self.vl_suffix is None:
            raise ValueError(
                "FitState cannot be instantiated until vl_ndata, tr_ndata and vl_suffix are filled"
            )
        self.training = training_info
        self.validation = validation_info
        self._parsed = False
        self._vl_chi2 = None
        self._tr_chi2 = None
        self._vl_dict = None
        self._tr_dict = None

    @property
    def vl_loss(self):
        """Return the total validation loss as it comes from the info dictionaries"""
        return self.validation.get("loss")

    @property
    def tr_loss(self):
        """Return the total validation loss as it comes from the info dictionaries"""
        return self.training.get("loss")

    def _parse_chi2(self):
        """
        Parses the chi2 from the losses according to the `tr_ndata` and
        `vl_ndata` dictionaries of {dataset: n_points}
        """
        if self._parsed:
            return
        if self.training is not None:
            self._tr_chi2, self._tr_dict = parse_losses(self.training, self.tr_ndata)
        if self.validation is not None:
            self._vl_chi2, self._vl_dict = parse_losses(
                self.validation, self.vl_ndata, suffix=self.vl_suffix
            )

    @property
    def tr_chi2(self):
        self._parse_chi2()
        return self._tr_chi2

    @property
    def vl_chi2(self):
        self._parse_chi2()
        return self._vl_chi2

    @property
    def all_tr_chi2(self):
        self._parse_chi2()
        return self._tr_dict

    @property
    def all_vl_chi2(self):
        self._parse_chi2()
        return self._vl_dict

    def all_tr_chi2_for_replica(self, r):
        """Return the tr chi2 per dataset for a given replica"""
        return {k: np.take(i, r) for k, i in self.all_tr_chi2.items()}

    def all_vl_chi2_for_replica(self, r):
        """Return the vl chi2 per dataset for a given replica"""
        return {k: np.take(i, r) for k, i in self.all_vl_chi2.items()}

    def total_partial_tr_chi2(self):
        """Return the tr chi2 summed over replicas per experiment"""
        return {k: np.sum(i) for k, i in self.all_tr_chi2.items()}

    def total_partial_vl_chi2(self):
        """Return the vl chi2 summed over replicas per experiment"""
        return {k: np.sum(i) for k, i in self.all_tr_chi2.items()}

    def total_tr_chi2(self):
        """Return the total tr chi2 summed over replicas"""
        return np.sum(self.tr_chi2)

    def total_vl_chi2(self):
        """Return the total vl chi2 summed over replicas"""
        return np.sum(self.vl_chi2)

    def __str__(self):
        return f"chi2: tr={self.tr_chi2} vl={self.vl_chi2}"


class ReplicaState:
    """Extra complication which eventually will be merged with someone else
    but it is here only for development."""

    def __init__(self, pdf_model):
        self._pdf_model = pdf_model
        self._weights = None
        self._best_vl_chi2 = INITIAL_CHI2

    @property
    def best_vl(self):
        return float(self._best_vl_chi2)

    def register_best(self, chi2, epoch):
        """Register a new best state and some metadata about it"""
        self._weights = self._pdf_model.get_weights()
        self._best_vl_chi2 = chi2

    def reload(self):
        """Reload the weights of the best state"""
        if self._weights:
            self._pdf_model.set_weights(self._weights)


class FitHistory:
    """
    Keeps a list of FitState items holding the full history of the fit.

    It also keeps track of the best epoch and the associated weights.

    Parameters
    ----------
        pdf_models: n3fit.backends.MetaModel
            list of PDF models being trained, used to saved the weights
    """

    def __init__(self, pdf_models, tr_ndata, vl_ndata):
        # Create a ReplicaState object for all models
        # which will hold the best chi2 and weights per replica
        self._replicas = []
        for pdf_model in pdf_models:
            self._replicas.append(ReplicaState(pdf_model))

        if vl_ndata is None:
            vl_ndata = tr_ndata
            vl_suffix = "loss"
        else:
            vl_suffix = "val_loss"
        # All instances of FitState should use these
        FitState.tr_ndata = tr_ndata
        FitState.vl_ndata = vl_ndata
        FitState.vl_suffix = vl_suffix

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

    def save_best_replica(self, i, epoch=None):
        """Save the state of replica ``i`` as a best fit so far.
        If an epoch is given, save the best as the given epoch, otherwise
        use the last one
        """
        if epoch is None:
            epoch = self.final_epoch
        loss = self.get_state(epoch).vl_loss[i]
        self._replicas[i].register_best(loss, epoch)

    def all_positivity_status(self):
        """Returns whether the positivity passed or not per replica"""
        return np.array([i.positivity_status for i in self._replicas])

    def all_best_vl_loss(self):
        """Returns the best validation loss for each replica"""
        return np.array([i.best_vl for i in self._replicas])

    def register(self, epoch, training_info, validation_info):
        """Save a new fitstate and updates the current final epoch

        Parameters
        ----------
            fitstate: FitState
                FitState object
                the fitstate of the object to save
            epoch: int
                the current epoch of the fit
        """
        # Save all the information in a fitstate object
        fitstate = FitState(training_info, validation_info)
        self.final_epoch = epoch
        self._history.append(fitstate)
        return fitstate

    def reload(self):
        """Reloads the best fit weights into the model if there are models to be reloaded
        Ensure that all replicas have stopped at this point.
        """
        for replica in self._replicas:
            replica.reload()


class Stopping:
    """
    Driver of the stopping algorithm

    Note, if the total number of points in the validation dictionary is None, it is assumed
    the validation_model actually corresponds to the training model.

    Parameters
    ----------
        validation_model: n3fit.backends.MetaModel
           the model with the validation mask applied
           (and compiled with the validation data and covmat)
        all_data_dict: dict
           list containg all dictionaries containing all information about
           the experiments/validation/regularizers/etc to be parsed by Stopping
        pdf_models: list(n3fit.backends.MetaModel)
           list of pdf_models being trained
        threshold_positivity: float
           maximum value allowed for the sum of all positivity losses
        total_epochs: int
           total number of epochs
        stopping_patience: int
           how many epochs to wait for the validation loss to improve
        dont_stop: bool
           dont care about early stopping
    """

    def __init__(
        self,
        validation_model,
        all_data_dicts,
        pdf_models,
        threshold_positivity=THRESHOLD_POS,
        total_epochs=0,
        stopping_patience=7000,
        threshold_chi2=10.0,
        dont_stop=False,
    ):
        # Save the validation object
        self._validation = validation_model

        # Create the History object
        tr_ndata, vl_ndata, pos_sets = parse_ndata(all_data_dicts)
        self._history = FitHistory(pdf_models, tr_ndata, vl_ndata)

        # And the positivity checker
        self._positivity = Positivity(threshold_positivity, pos_sets)

        # Initialize internal variables for the stopping
        self.n_replicas = len(pdf_models)
        self.threshold_chi2 = threshold_chi2
        self.stopping_degree = np.zeros(self.n_replicas, dtype=int)
        self.count = np.zeros(self.n_replicas, dtype=int)

        self.dont_stop = dont_stop
        self.stop_now = False
        self.stopping_patience = stopping_patience
        self.total_epochs = total_epochs

        self.stop_epochs = [None] * self.n_replicas
        self.best_epochs = [None] * self.n_replicas
        self.positivity_statusses = [POS_BAD] * self.n_replicas

    @property
    def vl_chi2(self):
        """Current validation chi2"""
        validation_info = self._validation.compute_losses()
        fitstate = FitState(None, validation_info)
        return fitstate.vl_chi2

    @property
    def e_best_chi2(self):
        """Epoch of the best chi2, if there is no best epoch, return last"""
        best_or_last_epochs = [
            best if best is not None else last
            for best, last in zip(self.best_epochs, self.stop_epochs)
        ]
        return best_or_last_epochs

    @property
    def stop_epoch(self):
        """Epoch in which the fit is stopped"""
        return self._history.final_epoch + 1

    @property
    def positivity_status(self):
        """Returns POS_PASS if positivity passes or veto if it doesn't
        for each replica"""
        return self._history.all_positivity_status()

    def evaluate_training(self, training_model):
        """Given the training model, evaluates the
        model and parses the chi2 of the training datasets

        Parameters
        ----------
            training_model: n3fit.backends.MetaModel
                an object implementing the evaluate function

        Returns
        -------
            tr_chi2: float
                chi2 of the given ``training_model``
        """
        training_info = training_model.compute_losses()
        fitstate = FitState(training_info, None)
        return fitstate.tr_chi2

    def monitor_chi2(self, training_info, epoch, print_stats=False):
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
            training_info: dict
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
        if np.isnan(training_info["loss"]):
            log.warning(" > NaN found, stopping activated")
            self.make_stop()
            return False

        # Step 2. Compute the validation metrics
        validation_info = self._validation.compute_losses()

        # Step 3. Register the current point in (the) history
        fitstate = self._history.register(epoch, training_info, validation_info)
        if print_stats:
            self.print_current_stats(epoch, fitstate)

        # Step 4. Check whether this is a better fit
        #         this means improving vl_chi2 and passing positivity
        # Don't start counting until the chi2 of the validation goes below a certain threshold
        # once we start counting, don't bother anymore
        passes = self.count | (fitstate.vl_chi2 < self.threshold_chi2)
        passes &= fitstate.vl_loss < self._history.all_best_vl_loss()
        # And the ones that pass positivity
        passes &= self._positivity(fitstate)

        self.stopping_degree += self.count

        # Step 5. loop over the valid indices to check whether the vl improved
        for i in np.where(passes)[0]:
            self.best_epochs[i] = epoch
            # By definition, if we have a ``best_epoch`` then positivity passed
            self.positivity_statusses[i] = POS_OK

            self._history.save_best_replica(i)
            self.stopping_degree[i] = 0
            self.count[i] = 1

        stop_replicas = self.count & (self.stopping_degree > self.stopping_patience)
        for i in np.where(stop_replicas)[0]:
            self.stop_epochs[i] = epoch
            self.count[i] = 0

        # By using the stopping degree we only stop when none of the replicas are improving anymore
        if min(self.stopping_degree) > self.stopping_patience:
            self.make_stop()
        return True

    def make_stop(self):
        """Convenience method to set the stop_now flag
        and reload the history to the point of the best model if any
        """
        self.stop_now = True
        self._history.reload()

    def print_current_stats(self, epoch, fitstate):
        """
        Prints ``fitstate`` training and validation chi2s
        """
        epoch_index = epoch + 1
        tr_chi2 = fitstate.total_tr_chi2()
        vl_chi2 = fitstate.total_vl_chi2()
        total_str = f"At epoch {epoch_index}/{self.total_epochs}, total chi2: {tr_chi2}\n"

        # The partial chi2 makes no sense for more than one replica at once:
        if self.n_replicas == 1:
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
        if self.dont_stop:
            return False
        else:
            return self.stop_now

    def chi2exps_json(self, replica=0, log_each=100):
        """
        Returns and apt-for-json dictionary with the status of the fit every `log_each` epochs

        Parameters
        ----------
            replica: int
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

        for i in range(log_each - 1, final_epoch + 1, log_each):
            fitstate = self._history.get_state(i)
            all_tr = fitstate.all_tr_chi2_for_replica(replica)
            all_vl = fitstate.all_vl_chi2_for_replica(replica)

            tmp = {exp: {"training": tr_chi2} for exp, tr_chi2 in all_tr.items()}
            for exp, vl_chi2 in all_vl.items():
                if exp not in tmp:
                    tmp[exp] = {"training": None}
                tmp[exp]["validation"] = vl_chi2

            json_dict[i + 1] = tmp
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

    def __call__(self, fitstate):
        """
        Checks whether a given FitState object
        passes the positivity requirement
        """
        return self.check_positivity(fitstate.validation)
