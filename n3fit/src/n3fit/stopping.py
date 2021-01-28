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
            all_tr_chi2: dict
                Chi2 for all training datasets computed before the update of the weights
            all_vl_chi2: dict
                Chi2 for all validation datasets computed after the update of the weights
            info: dict
                Full return state of the Validation model
    """

    def __init__(self, all_tr_chi2, all_vl_chi2, info):
        self.all_tr_chi2 = all_tr_chi2
        self.all_vl_chi2 = all_vl_chi2
        self.info = info
        # These two variables are only filled for specific points
        # in order to save precious memory, and only when we are
        # saving the fit history each X number of epoch
        self.weights = None
        self.best_epoch = 0

    @property
    def vl_chi2(self):
        """ Returns the total validation chi2 """
        return self.all_vl_chi2["total"]

    @property
    def tr_chi2(self):
        """ Returns the total training chi2 """
        return self.all_tr_chi2["total"]

    def vl_chi2_for_replica(self, i):
        """ Returns the validation_chi2 for a given replica """
        return self.all_vl_chi2["total"][i]

    def __str__(self):
        return f"chi2: tr={self.tr_chi2} vl={self.vl_chi2}"


class ReplicaBest:
    """ Extra complication which eventually will be merged with someone else
    but it is here only for development."""

    def __init__(self, pdf_model):
        self._pdf_model = pdf_model
        self._weights = None
        self._best_epoch = None
        self._best_vl_chi2 = INITIAL_CHI2

    def positivity_pass(self):
        """ By definition, if we have a ``best_epoch`` then positivity passed """
        if self._best_epoch is None:
            return False
        else:
            return True

    @property
    def best_epoch(self):
        return self._best_epoch

    @property
    def best_vl(self):
        return float(self._best_vl_chi2)

    @property
    def positivity_status(self):
        if self.positivity_pass():
            return POS_OK
        else:
            return POS_BAD

    def register_best(self, chi2, epoch):
        self._weights = self._pdf_model.get_weights()
        self._best_epoch = epoch
        self._best_vl_chi2 = chi2

    def reload(self):
        if self._weights:
            self._pdf_model.set_weights(self._weights)


class FitHistory:
    """
        Keeps a list of FitState items holding the full history of the fit.

        It also keeps track of the best epoch and the associated weights.

        Can be iterated when there are snapshots of the fit being saved.
        When iterated it will rewind the fit to each of the point in history
        that have been saved.

        Parameters
        ----------
            pdf_models: n3fit.backends.MetaModel
                list of PDF models being trained, used to saved the weights

            save_weights_each: int
                if given, it will save a snapshot of the fit every  `save_weights_each` epochs
    """
    def __init__(self, pdf_models, save_weights_each=None):
        # Save a list of status per replica
        self._replicas = []
        for pdf_model in pdf_models:
            self._replicas.append(ReplicaBest(pdf_model))
        # Save a list of status for the entire fit
        self._history = []

        self._save_weights_each = save_weights_each
        # Initialize variables for the history
        self._weights = None
        self._best_epoch = None
        self.final_epoch = None
        # Initialize variables for the snapshots
        self.reloadable_history = []

    @property
    def best_epoch(self):
        """ Epoch of the best fit """
        return self._best_epoch

    @best_epoch.setter
    def best_epoch(self, epoch):
        """ Saves the current weight """
        self._weights = self._pdf_model.get_weights()
        self._best_epoch = epoch

    def get_state(self, epoch):
        """ Get the FitState of the system for a given epoch """
        return self._history[epoch]

    def best_state(self):
        """ Return the FitState object corresponding to the best fit """
        if self.best_epoch is None:
            return None
        else:
            index = self.best_epoch
            best_state = self._history[index]
            return best_state

    def best_vl(self):
        """ Returns the chi2 of the best fit
        if there was no best fit returns `INITIAL_CHI2`
        if there was a problem, returns `TERRIBLE_CHI2` """
        if not self._weights:
            return TERRIBLE_CHI2
        best_state = self.best_state()
        if best_state:
            return best_state.vl_chi2
        else:
            return INITIAL_CHI2

    def best_tr(self):
        """ Returns the training chi2 of the best fit
        if there are no best fit, returns the last one """
        best_state = self.best_state()
        if best_state:
            return best_state.tr_chi2
        else:
            return self._history[self.final_epoch].tr_chi2

    def save_best_replica(self, i, epoch = None):
        """ Save the state of replica ``i`` as a best fit so far.
        If an epoch is given, save the best as the given epoch, otherwise
        use the last one
        """
        if epoch is None:
            epoch = self.final_epoch
        chi2 = self._history[epoch].vl_chi2[i]
        self._replicas[i].register_best(chi2, epoch)

    def all_best_vl(self):
        """ Returns the best validation chi2 for each replica """
        return [i.best_vl for i in self._replicas]

    def register(self, fitstate, epoch):
        """ Save a new fitstate and updates the current final epoch
        Every `save_weights_each` (if set) saves a snapshot of the current best fit into
        the fitstate

        Parameters
        ----------
            fitstate: FitState
                FitState object 
                the fitstate of the object to save
            `epoch`
                the current epoch of the fit
        """
        self.final_epoch = epoch
        self._history.append(fitstate)
        if self._save_weights_each:
            save_here = (epoch + 1) % self._save_weights_each
            # TODO this must be done differently now
#             if save_here == 0:
#                 fitstate.register_weigths(self._weights, self.best_epoch)
#                 self.reloadable_history.append(fitstate)

    def reload(self, weights=None):
        """ Reloads the best fit weights into the model
        if there are models to be reloaded
        A set of weights can be enforced as an optional argument
        """
        # TODO the weights part is tricky
        if weights is None:
            for replica in self._replicas:
                replica.reload()
#         if weights:
#             self._pdf_model.set_weights(weights)

    def rewind(self, step):
        """ Rewind the FitHistory object to the step `step` in the fit
        """
        fitstate = self.reloadable_history[step]
        historic_weights = fitstate.weights
        self.reload(weights=historic_weights)
        self.best_epoch = fitstate.best_epoch
        self.final_epoch = (step + 1) * self._save_weights_each - 1
        # -1 because we are saving the epochs starting at 0


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
            save_weights_each: int
                every how many epochs to save a snapshot of the fit
    """

    def __init__(
        self,
        validation_model,
        all_data_dicts,
        pdf_models,
        threshold_positivity=1e-6,
        total_epochs=0,
        stopping_patience=7000,
        threshold_chi2=10.0,
        dont_stop=False,
        save_weights_each=None,
    ):
        # Parse the training, validation and positivity sets from all the input dictionaries
        self._tr_ndata, vl_ndata, pos_sets = parse_ndata(all_data_dicts)

        # Create the Validation, Positivity and History objects
        if vl_ndata is None:
            self.validation = Validation(validation_model, self._tr_ndata, no_validation=True)
        else:
            self.validation = Validation(validation_model, vl_ndata)
        self.positivity = Positivity(threshold_positivity, pos_sets)

        self.history = FitHistory(pdf_models)

        # Initialize internal variables for the stopping
        self.n_replicas = len(pdf_models)
        self.threshold_chi2 = threshold_chi2
        self.dont_stop = dont_stop
        self.stop_now = False
        self.stopping_patience = stopping_patience
        self.stopping_degree = 0
        self.count = 0
        self.total_epochs = total_epochs
        self.replica_iterator = None

    @property
    def vl_chi2(self):
        """ Validation chi2 """
        return self.history.best_vl()

    @property
    def tr_chi2(self):
        """ Training chi2 """
        return self.history.best_tr()

    @property
    def e_best_chi2(self):
        """ Epoch of the best chi2, if there is no best epoch
        return the last epoch"""
        be = self.history.best_epoch
        if be is None:
            return self.stop_epoch
        return be


    @property
    def stop_epoch(self):
        """ Epoch in which the fit is stopped """
        return self.history.final_epoch + 1

    def evaluate_training(self, training_model):
        """ Given the training model, evaluates the
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
        tr_chi2, _ = parse_losses(training_info, self._tr_ndata)
        return tr_chi2

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
        # Step 1. Preprocess the event, count it towards the stopping degree
        #         parse the training information and check whether it is a good point
        tr_chi2, all_tr = parse_losses(training_info, self._tr_ndata)
        
        if np.isnan(tr_chi2):
            log.warning(" > NaN found, stopping activated")
            self.stop_now = True
            # If we had a good model at any point, reload
            self.history.reload() # TODO
            return False

        self.stopping_degree += self.count

        # Step 2. Check the validation loss at this point
        # each loss is an array of the loss per replica
        vl_chi2, all_vl = self.validation.loss()

        # Step 3. Store information about the run and print stats if asked
        fitstate = FitState(all_tr, all_vl, self.validation.state)
        self.history.register(fitstate, epoch)

        if print_stats:
            self.print_current_stats(epoch, fitstate)

        # Step 4. Check whether this is a better fit
        #         this means improving vl_chi2 and passing positivity

        # Get the values that pass validation
        passes_val = vl_chi2 < self.threshold_chi2
        passes_val &= vl_chi2 < self.history.all_best_vl()
        # And the ones that pass positivity
        passes_pos = self.positivity(fitstate)

        # Now loop over the valid indices to check whether the vl improved
        # TODO: check whether this loop is hurting at all performance (shouldnt???)
        for i in np.where(passes_val & passes_pos)[0]:
            self.history.save_best_replica(i)

            # There is no stopping for now
        return True
        
        # TODO for now we force all fits to get to the end


        # For each of them, check whether the validation chi2 is better or not
        if self.positivity(fitstate) and vl_chi2 < self.threshold_chi2:
            if vl_chi2 < self.history.best_vl():
                # Set the new best
                self.history.best_epoch = epoch
                # Reset stopping info
                self.stopping_degree = 0
                # Initialize the counter
                self.count = 1

        # If your patience has ended, prepare for stop
        if self.stopping_degree > self.stopping_patience:
            self.make_stop()
        return True

    def make_stop(self):
        """ Convenience method to set the stop_now flag
        and reload the history to the point of the best model if any
        """
        self.stop_now = True
        self.history.reload()

    def print_current_stats(self, epoch, fitstate):
        """
            Prints the last validation and training loss saved
        """
        epoch_index = epoch + 1
        tr_loss = fitstate.tr_chi2
        vl_loss = fitstate.vl_chi2
        total_str = f"At epoch {epoch_index}/{self.total_epochs}, total loss: {tr_loss}\n"

        partials = []
        for experiment in self._tr_ndata:
            chi2 = fitstate.all_tr_chi2[experiment]
            partials.append(f"{experiment}: {chi2:.3f}")
        total_str += ", ".join(partials)

        total_str += f"\nValidation loss at this point: {vl_loss}"
        log.info(total_str)

    def stop_here(self):
        """ Returns the stopping status
        If `dont_stop` is set returns always False (i.e., never stop)
        """
        if self.dont_stop:
            return False
        else:
            return self.stop_now

    def positivity_pass(self):
        """ Checks whether the positivity loss is below the requested threshold
        If there is no best state, the positivity (obv) cannot pass
        """
        best_state = self.history.best_state()
        if best_state is not None and self.positivity(best_state):
            return True
        else:
            return False

    def positivity_status(self):
        """ Checks whether the positivity loss is below the requested threshold
        If there is no best state, the positivity (obv) cannot pass
        """
        if self.positivity_pass():
            return POS_OK
        else:
            return POS_BAD

    def get_next_replica(self):
        """ Return the next ReplicaBest object"""
        if self.replica_iterator is None:
            self.replica_iterator = iter(self.history._replicas)
            self.ii = -1
        self.ii += 1
        return self.ii, next(self.replica_iterator)

    def chi2exps_str(self, log_each=100):
        """
        Returns a list of log-string with the status of the fit
        every `log_each` epochs

        Parameters
        ----------
            `log_each`
                every how many epochs to print the log

        Returns
        -------
            `file_list`
                a list of string to be printed as `chi2exps.log`
        """
        final_epoch = self.history.final_epoch
        file_list = []
        for i in range(log_each - 1, final_epoch + 1, log_each):
            fitstate = self.history.get_state(i)
            all_tr = fitstate.all_tr_chi2
            all_vl = fitstate.all_vl_chi2
            # Here it is assumed the validation exp set is always a subset of the training exp set
            data_list = []
            for exp in self._tr_ndata:
                tr_loss = all_tr[exp]
                vl_loss = all_vl.get(exp, 0.0)
                data_str = f"{exp}: {tr_loss} {vl_loss}"
                data_list.append(data_str)
            data = "\n".join(data_list)
            epoch_index = i + 1
            total_tr_loss = fitstate.tr_chi2
            total_vl_loss = fitstate.vl_chi2
            strout = f"""
Epoch: {epoch_index}
{data}
Total: training = {total_tr_loss} validation = {total_vl_loss}
"""
            file_list.append(strout)
        return file_list


class Validation:
    """
        Controls the NNPDF cross-validation algorithm

        The cross-validation refers to the validation loss of the points of the dataset
        not used in the fitting.
        In general for any points considered here there will accompanying points from the
        same dataset being included in the fitting.

        Parameters
        ----------
            model: n3fit.backends.MetaModel
                the model with the validation mask applied
                (and compiled with the validation data and covmat)
    """

    def __init__(self, model, ndata_dict, verbose=False, no_validation=False):
        self.model = model
        self.state = None
        self.verbose = verbose
        self.ndata_dict = ndata_dict
        self.n_val_exp = len(ndata_dict)
        if no_validation:
            self.suffix = "loss"
        else:
            self.suffix = "val_loss"

    def _compute_validation_loss(self):
        """
        Evaluates the validation model and returns a tuple (`total_loss`, `vl_dict`)
        with the information for the validation loss by experimenet normalized to the
        number of points of each experiment

        Returns
        -------
            total_loss: float
                total vale for the validation loss
            vl_dict: dict
                dictionary containing a map of experiment names and their loss per replica
        """
        loss_dict = self.model.compute_losses()
        self.state = loss_dict
        return parse_losses(loss_dict, self.ndata_dict, suffix=self.suffix)

    def loss(self):
        """
        Returns a tuple with the validation loss and a
        dictionary for the validation loss per experiment
        """
        return self._compute_validation_loss()


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

        Parameters
        ----------
            history_object: dict
                dictionary of entries in the form  {'name': loss}, output of a MetaModel .fit()
        """
        positivity_pass = True
        for key in self.positivity_sets:
            key_loss = f"{key}_loss"
            positivity_pass &= history_object[key_loss] < self.threshold
            if not all(positivity_pass):
                return np.array(positivity_pass)
        return np.array(positivity_pass)

    def __call__(self, fitstate):
        """
            Checks whether a given FitState object
            passes the positivity requirement
        """
        return self.check_positivity(fitstate.info)
