"""
    Module containing the classes related to the stopping alogirthm

    In this module there are three Classes:
    - FitHistory: this class contains the information necessary
            in order to reset the state of the fit to the point
            in which the history was saved.
    - Stopping: this class monitors the chi2 of the validation
            and training sets and decides when to stop
    - Positivity: Decides whether a given point fullfills the positivity conditions
"""

# TODO
#   - Does it make sense to use Keras callbacks or we want to actually avoid it?
#           Extra: do we want to implement the callback as part of the .fit function?

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

    # Return:
        - `tr_ndata`: dictionary of {'exp' : ndata}
        - `vl_ndata`: dictionary of {'exp' : ndata}
        - `pos_set`: list of the names of the positivity sets

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
        if dictionary.get("positivity"):
            pos_set.append(exp_name)
    if not vl_ndata_dict:
        vl_ndata_dict = tr_ndata_dict
    return tr_ndata_dict, vl_ndata_dict, pos_set

class FitState:
    """
        Holds the state of the fit for a given point in history
    """
    def __init__(self, all_tr_chi2, all_vl_chi2, training_info):
        self.all_tr_chi2 = all_tr_chi2
        self.all_vl_chi2 = all_vl_chi2
        self.training_info = training_info
        # These two variables are only filled for specific points
        # in order to save precious memory, and only when we are
        # saving the fit history each X number of epoch
        self.weights = None
        self.best_epoch = 0

    @property
    def vl_chi2(self):
        return self.all_vl_chi2["total"]

    @property
    def tr_chi2(self):
        return self.all_tr_chi2["total"]

    def __str__(self):
        return f"chi2: tr={self.tr_chi2} vl={self.vl_chi2}"

    def save_history(self, weights, best_epoch):
        self.weights = weights
        self.best_epoch = best_epoch

class FitHistory:
    """
        List of FitState items.
        It also keeps track of which epoch produced the best fit and the associated weights
    """
    def __init__(self, validation_model, save_each = None):
        self.final_epoch = 0
        self.best_epoch = -1
        self.weights = None
        self.history = []
        self.reloadable_history = []
        self.validation_model = validation_model
        self.save_each = save_each
        self.terrible = False

    @property
    def best_state(self):
        if self.best_epoch < 0:
            return None
        else:
            index = self.best_epoch
            best_state = self.history[index]
            return best_state

    def best_vl(self):
        if self.terrible:
            return TERRIBLE_CHI2
        if self.best_state:
            return self.best_state.vl_chi2
        else:
            return INITIAL_CHI2

    def best_tr(self):
        if self.best_state:
            return self.best_state.tr_chi2
        else:
            return self.history[-1].tr_chi2

    def save(self, fitstate, epoch):
        self.final_epoch = epoch
        self.history.append(fitstate)
        if self.save_each:
            save_here = (epoch + 1) % self.save_each
            if save_here == 0:
                fitstate.save_history(self.weights, self.best_epoch)
                self.reloadable_history.append(fitstate)

    def new_best(self, epoch):
        self.best_epoch = epoch
        self.weights = self.validation_model.weights

    def reload(self):
        if self.weights:
            self.validation_model.weights = self.weights
        else:
            # If there was no model at this point, this was
            # a terrible run, mark it as such
            self.terrible = True


    def __iter__(self):
        for i, fitstate in enumerate(self.reloadable_history):
            log.info("Reloading step %d", i)
            self.weights = fitstate.weights
            self.best_epoch = fitstate.best_epoch
            self.final_epoch = (i+1)*self.save_each
            # TODO reload the log also so we have a partial log of chi2
            self.reload()
            yield i




class Stopping:
    """
        Driver of the stopping algorithm

        Note, if the total number of points in the validation dictionary is None, it is assumed
        the validation_model actually corresponds to the training model.

        # Arguments:
            - `validation_model`: the model with the validation mask applied
                                  (and compiled with the validation data and covmat)
            - `all_data_dict`: list containg all dictionaries containing all information about
                              the experiments/validation/regularizers/etc to be parsed by Stopping
            - `threshold_positivity`: maximum value allowed for the sum of all positivity losses
            - `total_epochs`: total number of epochs
            - `stopping_patience`: how many epochs to wait for the validation loss to improve
            - `dont_stop`: dont care about early stopping

    """

    def __init__(
        self,
        validation_model,
        all_data_dicts,
        threshold_positivity=1e-6,
        total_epochs=0,
        stopping_patience=7000,
        dont_stop=False,
        save_each=None,
    ):
        # Parse from the 
        self.ndata_tr_dict, vl_ndata, pos_sets = parse_ndata(all_data_dicts)

        # Save the validation and positivity objects
        self.validation = Validation(validation_model, vl_ndata)
        self.positivity = Positivity(threshold_positivity, pos_sets)

        # Initialize internal (useful) variables
        self.save_history = False
        self.stop_now = False
        self.nan_found = False
        self.dont_stop = dont_stop
        self.stopping_patience = stopping_patience
        self.stopping_degree = 0
        self.count = 0
        self.partial_log = -1

        # Initialize the parameters of the fit
        self.total_epochs = total_epochs
        self.epoch_of_the_stop = total_epochs
        self.training_chi2 = INITIAL_CHI2
        self.best_chi2 = INITIAL_CHI2
        self.history = FitHistory(self.validation)

        # Initialize stats
        self.reset_stats()

    # In order to get the validation and training loss we return the loss of the best epoch
    # the reason we do this is the HistoryMaker, the easiest and cleanest way of rewinding history
    # is to go back to the epoch in which the weights for that point in history were saved
    @property
    def vl_loss(self):
        """ Validation loss """
        return self.history.best_vl()

    @property
    def tr_loss(self):
        """ Training loss """
        return self.history.best_tr()

    @property
    def e_best_chi2(self):
        """ Epoch of the best chi2 """
        return self.history.best_epoch

    def reset_stats(self):
        """ Reset statistical information """
        self.training_losses = []
        self.validation_losses = []
        self.file_list = []
        self.w_best_chi2 = None

    def history_loop(self):
        """
        Loop over the history object
        If history is active will return an iterable object ranging from 0 to n_steps saved
        Each iteration of the loop reloads the corresponding point in history (weights and chi2s)
        """
        if self.save_history:
            return self.history
        else:
            return []

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

        # Arguments:
            - `training_info`: the output of a .fit() run
            - `epoch`: the index of the epoch

        # Returns:
            - `pass_ok`: true/false according to the status of the run
        """
        # Step 1. Preprocess the event, count it towards the stopping degree
        #         parse the training information and check whether it is a good point
        tr_chi2, all_tr = self._parse_training(training_info)

        if np.isnan(tr_chi2):
            log.warning(" > NaN found, stopping activated")
            self.stop_now = True
            self.nan_found = True
            self.epoch_of_the_stop = epoch
            # If we had a good model at any point, reload
            self.history.reload()
            return False

        self.stopping_degree += self.count

        # Step 2. Check the validation loss at this point
        vl_chi2, all_vl = self.validation.loss()

        # Step 3. Store information about the run
        fitstate = FitState(all_tr, all_vl, training_info)
        self.history.save(fitstate, epoch)

        # Step 4. Check whether this is a better fit
        #         this means improving vl_chi2 and passing positivity
        if self.positivity(fitstate):
            if vl_chi2 < self.history.best_vl():
                # Set the new best
                self.history.new_best(epoch)
                # Save stopping info
                self.stopping_degree = 0
                # Initialize the counter
                self.count = 1

        # If we are asked, print stats and save logs
        if print_stats:
            self.save_chi2_log(epoch, fitstate)
            self.print_current_stats(epoch, fitstate)

        # If your patience has ended, prepare for stop
        if self.stopping_degree > self.stopping_patience:
            self.stop_now = True
            self.epoch_of_the_stop = epoch
            self.history.reload()
        return True

    def print_current_stats(self, epoch, fitstate):
        """
            Prints the last validation and training loss saved
        """
        epoch_index = epoch + 1
        tr_loss = fitstate.tr_chi2
        vl_loss = fitstate.vl_chi2
        total_str = f"At epoch {epoch_index}/{self.total_epochs}, total loss: {tr_loss}\n"

        partials = []
        for experiment in self.ndata_tr_dict:
            chi2 = fitstate.all_tr_chi2[experiment]
            partials.append(f"{experiment}: {chi2:.3f}")
        total_str += ", ".join(partials)

        total_str += f"\nValidation loss at this point: {vl_loss}"
        log.info(total_str)

    def save_chi2_log(self, epoch, fitstate):
        # TODO remove this function as it is 100% useless
        """
        For each call save a string in a list with the information of the real current chi2
        (i.e., pays no attention to validation or positivity)
        This is used to generate a log file so that the evolution of the chi2 per experiment
        for the fit can be observed.

        # Arguments:
            - `epoch`: the current epoch
            - `all_tr`: dictionary of experiment names and their corresponding chi2 (training)
            - `all_vl`: dictionary of experiment names and their corresponding chi2 (validation)
        """
        all_tr = fitstate.all_tr_chi2
        all_vl = fitstate.all_vl_chi2
        # Note: here it is assumed the validation exp set is always a subset of the training exp set
        data_list = []
        for exp, tr_loss in all_tr.items():
            vl_loss = all_vl.get(exp, 0.0)
            data_str = f"{exp}: {tr_loss} {vl_loss}"
            data_list.append(data_str)
        data = "\n".join(data_list)
        epoch_index = epoch + 1
        total_tr_loss = fitstate.tr_chi2
        total_vl_loss = fitstate.vl_chi2
        strout = f"""
Epoch: {epoch_index}
{data}
Total: training = {total_tr_loss} validation = {total_vl_loss}
"""
        self.file_list.append(strout)


    def _parse_training(self, training_info):
        """
        Receives an object containg the training chi2.
        Usually a history object, but it can come in the form of a dictionary.

        It loops over the dictionary and uses the npoints_data dictionary to
        normalize the chi2 and return backs a tuple (`total`, `tr_chi2`)

        # Arguments:
            - `training_info`: history object

        # Returns:
            - `total` : total value for the training loss
            - `tr_chi2`: dictionary of {'expname' : loss }
        """
        try:
            hobj = training_info.history
        except AttributeError:  # So it works whether we pass the out our the out.history
            hobj = training_info

        # In the general case epochs = 1.
        # In case that we are doing more than 1 epoch, take the average to smooth out
        # fluctuations.
        # This value is only used for printing output purposes so should not have any significance
        tr_chi2 = {}
        total_points = 0
        total_loss = 0
        for exp_name, npoints in self.ndata_tr_dict.items():
            loss = np.mean(hobj[exp_name + "_loss"])
            tr_chi2[exp_name] = loss / npoints
            total_points += npoints
            total_loss += loss

        # By taking the loss from the history object we would be saving the total loss
        # including positivity sets and (if added/enabled) regularizsers
        # instead we want to restrict ourselves to the loss coming from experiments
        # total_loss = np.mean(hobj["loss"]) / total_points
        total_loss /= total_points
        tr_chi2["total"] = total_loss
        return total_loss, tr_chi2

    def stop_here(self):
        """ Returns the stopping status unless dont_stop is true, then returns False """
        if self.dont_stop:
            return False
        else:
            return self.stop_now

    def positivity_pass(self):
        """ Checks whether the positivity loss is below the requested threshold """
        if self.positivity(self.history.best_state):
            return POS_OK
        else:
            return POS_BAD

    def chi2exps_str(self):
        """
        Returns the list of log-strings.
        If history has been reloaded, there will be a variable `partial_log` and
        only the log up to that point will be returned
        """
        return self.file_list[: self.partial_log]


class Validation:
    """
        Controls the NNPDF cross-validation algorithm

        The cross-validation refers to the validation loss of the points of the dataset
        not used in the fitting.
        In general for any points considered here there will accompanying points from the
        same dataset being included in the fitting.

        # Arguments:
            - `model`: the model with the validation mask applied
                       (and compiled with the validation data and covmat)
    """

    def __init__(self, model, ndata_dict, verbose=False):
        self.model = model
        self.verbose = verbose
        self.ndata_dict = ndata_dict
        # If there are extra losses they will appear at the end of the list, so we want to restrict
        # ourselves to the chi2, which means we want to go up to the number of exp. with validation
        self.n_val_exp = len(ndata_dict)

    def _compute_validation_loss(self):
        """
        Evaluates the validation model and returns a tuple (`total_loss`, `vl_dict`)
        with the information for the validation loss by experimenet normalized to the
        number of points of each experiment

        # Returns:
            - `total_loss`: total vale for the validation loss
            - `vl_dict`: dictionary containing a map of experiment names and loss
        """
        # The variable vl_list is a list of all losses of the model, where the first element
        # is sum of all other elements
        loss_list = self.model.evaluate(verbose=self.verbose)

        # This loop relies on the information that comes through the input dict to be accurate
        # because since (at the moment) the list that evaluate returns has no names, we need to
        # assume they come in the correct order (same order as the traiing losses)
        vl_dict = {}
        total_points = 0
        total_loss = 0
        for loss, (exp_name, npoints) in zip(
            loss_list[1 : 1 + self.n_val_exp], self.ndata_dict.items()
        ):
            vl_dict[exp_name] = loss / npoints
            total_loss += loss
            total_points += npoints

        total_loss /= total_points
        vl_dict["total"] = total_loss

        return total_loss, vl_dict

    @property
    def weights(self):
        """ Returns the weights of the validation model """
        return self.model.get_weights()

    @weights.setter
    def weights(self, weights):
        """  Sets the weights on the validation model """
        self.model.set_weights(weights)

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

        # Arguments:
            - `threshold_positivity`: maximum value allowed for the sum of all positivity losses
    """

    def __init__(self, threshold, positivity_sets):
        self.threshold = threshold
        self.positivity_ever_passed = False
        self.positivity_sets = positivity_sets

    def check_positivity(self, history_object):
        """
            This function receives a history object and look for entries
            with the keyname: pos_key_something


            # Arguments:
                - `history_object`: a dictionary of entries in the form
                    {'name': loss}, output of a MetaModel .fit()
                - `pos_key`: `key that searchs for the positivity`
        """
        positivity_loss = 0.0
        for key in self.positivity_sets:
            key_loss = f"{key}_loss"
            positivity_loss = history_object[key_loss][-1]
        if positivity_loss > self.threshold:
            return False
        else:
            self.positivity_ever_passed = True
            return True

    def __call__(self, fitstate):
        """
            Checks whether a given FitState object
            passes the positivity requirement
        """
        return self.check_positivity(fitstate.training_info)
