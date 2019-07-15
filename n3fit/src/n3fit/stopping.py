"""
    Module containing the classes related to the stopping alogirhtm
"""

# TODO
#   - Save the chi2 logs in a better way
#   - Does it make sense to use Keras callbacks or we want to actually avoid it?
#           Extra: do we want to implement the callback as part of the .fit function?
#   - Save the animations

import numpy as np

import logging

log = logging.getLogger(__name__)

# Put a very big number here so that we for sure discard this run
# AND we have a clear marker that something went wrong, not just a bad fit
TERRIBLE_CHI2 = 1e10
INITIAL_CHI2 = 1e9

# Pass/veto keys
POS_OK = "POS_PASS"
POS_BAD = "POS_VETO"

def parse_ndata(ndata_dict):
    """ 
    Parses the ndata dictionary into training a validation
    dictionaries

    Ignores entires for "total" and positivity sets.

    # Return:
        - `tr_ndata`: dictionary of {'exp' : ndata}
        - `vl_ndata`: dictionary of {'exp' : ndata}

    Note: if the total number of validation points remain at 0
    it is understood there is no validation and the output vl_data will also point to tr_ndata
    """
    tr_ndata = {}
    vl_ndata = {}
    val_points = 0
    for key, npoints in ndata_dict.items():
        if key.lower().startswith(("pos", "total")):
            continue
        if key.lower().endswith("_vl"):
            vl_ndata[key] = npoints
            val_points += npoints
        else:
            tr_ndata[key] = npoints
    if npoints == 0:
        vl_ndata = tr_ndata
    return tr_ndata, vl_ndata

class Stopping:
    """
        Driver of the stopping algorithm

        Note, if the total number of points in the validation dictionary is None, it is assumed
        the validation_model actually corresponds to the training model.

        # Arguments:
            - `validation_model`: the model with the validation mask applied (and compiled with the validation data and covmat)
            - `ndata_dict`: a dictionary with the number of points of each experiment in the format:
                `exp`: npoints in the training model data
                `exp_vl`: npoints in the validation model data
                This dictionary is used to print the loss normalized to the number of points
            - `threshold_positivity`: maximum value allowed for the sum of all positivity losses
            - `total_epochs`: total number of epochs
            - `stopping_patience`: how many epochs to wait for the validation loss to improve
            - `dont_stop`: dont care about early stopping

    """

    def __init__(
        self,
        validation_model,
        ndata_dict,
        threshold_positivity=1e-6,
        total_epochs=0,
        stopping_patience=7000,
        dont_stop=False,
    ):

        # Parse the data dictionary
        self.ndata_tr_dict, vl_ndata = parse_ndata(ndata_dict)

        # Save the validation and positivity objects
        self.validation = Validation(validation_model, vl_ndata)
        self.positivity = Positivity(threshold_positivity)

        # Initialize internal (useful) variables
        self.stop_now = False
        self.nan_found = False
        self.dont_stop = dont_stop
        self.stopping_patience = stopping_patience
        self.stopping_degree = 0
        self.count = 0

        # Initialize the parameters of the fit
        self.total_epochs = total_epochs
        self.epoch_of_the_stop = total_epochs
        self.best_chi2 = INITIAL_CHI2
        self.training_chi2 = INITIAL_CHI2

        # Initialize stats
        self.reset_stats()

    @property
    def loss(self):
        return self.best_chi2
    @property
    def tr_loss(self):
        return self.training_chi2

    def reset_stats(self):
        self.training_losses = []
        self.validation_losses = []
        self.w_best_chi2 = None
        self.e_best_chi2 = 0

    def save_stats(self, tr_loss, vl_loss):
        """ 
        Saves the training and validation losses each epoch

        # Argument:
         - `tr_loss`: training loss
         - `vl_loss`: validation loss
        """
        self.training_losses.append(tr_loss)
        self.validation_losses.append(vl_loss)

    def monitor_chi2(self, training_info, epoch, print_stats = False):
        """
        Function to be called at the end of every epoch.
        Stores the total chi2 of the training set as well as the
        total chi2 of the validation set.
        If the training chi2 is below a certain threshold,
        stores the state of the model which gave the minimum chi2
        as well as the epoch in which occurred
        If the epoch is a multiple of save_all_each then we also save the per-exp chi2

        Returns True if the run seems ok and False if a NaN is found

        # ArgumentS:
            - `training_info`: the output of a .fit() run
            - `epoch`: the index of the epoch

        # Returns:
            - `pass_ok`: true/false according to the status of the run
        """
        # Step 0. Count the event, worst case scenario we have to reset to 0
        self.stopping_degree += self.count

        # Step 1. Check whether the run was ok
        tr_chi2, all_tr = self._parse_training(training_info)

        if np.isnan(tr_chi2):
            log.warning(" > NaN found, stopping activated")
            self.stop_now = True
            self.nan_found = True
            self.epoch_of_the_stop = epoch
            # If we had a good model at any point, reload
            if self.w_best_chi2:
                self.reload_model()
            else:  # else, make sure the loss is marked as terrible
                self.best_chi2 = TERRIBLE_CHI2
            return False

        # Step 2. Read the validation loss (if there is no validation loss it will correspond to the training)
        vl_chi2, _ = self.validation.loss

        # Step 3. Store all information about the run
        self.save_stats(tr_chi2, vl_chi2)

        # Step 4. Check whether positivity passes, if it doesn't ignore the point
        if self.positivity.check_positivity(training_info):
            if vl_chi2 < self.best_chi2:
                # we found a better validation loss! Reset All
                self.best_chi2 = vl_chi2
                # Save the training chi2 for this point
                self.training_chi2 = tr_chi2
                # Save extra info
                self.stopping_degree = 0
                self.e_best_chi2 = epoch
                self.best_history = training_info
                self.save_model()
                # Now activate the counter
                self.count = 1

        # If it makes sense, print states
        if print_stats:
            self.print_current_stats(all_tr)

        # If your patience has ended, prepare for stop
        if self.stopping_degree > self.stopping_patience:
            self.stop_now = True
            self.epoch_of_the_stop = epoch
            self.reload_model(force=True)

    def save_model(self):
        """
            Saves the weights of the NN for the validation model
        """
        self.w_best_chi2 = self.validation.weights

    def reload_model(self, force=False):
        """
            Reload the weights of the NN onto the validation model
        """
        if self.w_best_chi2 and (self.not_reload or force):
            self.validation.weights = self.w_best_chi2

    def print_current_stats(self, all_tr):
        """
            Prints the last validation and training loss saved
        """
        epoch = len(self.training_losses)
        loss = self.training_losses[-1]
        vl_loss = self.validation_losses[-1]
        total_str = "\nAt epoch {0}/{1}, total loss: {2}".format(
            epoch, self.total_epochs, loss
        )

        partials = ""
        for experiment, chi2 in all_tr.items():
            partials += "{0}: {1:.3f}".format(experiment, chi2)

        print(total_str)
        if partials:
            print(partials)

        print(" > Validation loss at this point: {0}".format(vl_loss))

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
        for exp_name, npoints in self.ndata_tr_dict.items():
            chi2 = np.mean(hobj[exp_name + "_loss"]) / npoints
            tr_chi2[exp_name] = chi2
            total_points = npoints

        total = np.mean(hobj["loss"]) / total_points
        return total, tr_chi2

    def stop_here(self):
        if self.dont_stop:
            return False
        else:
            return self.stop_now

    def positivity_pass(self):
        return self.positivity()

    def chi2exps_str(self):
        file_list = []
        for epoch in range(0, len(self.training_losses), 100):
            strout = "\nEpoch: {0}".format(epoch)
            strout += "\nTotal: training = {0} validation = {1}\n".format(self.training_losses[epoch], self.validation_losses[epoch])
            file_list.append(strout)
        return file_list


class Validation:
    """
        Controls the NNPDF cross-validation algorithm

        The cross-validation refers to the validation loss of the points of the dataset
        not used in the fitting.
        In general for any points considered here there will accompanying points from the 
        same dataset being included in the fitting.

        # Arguments:
            - `model`: the model with the validation mask applied (and compiled with the validation data and covmat)
    """

    def __init__(self, model, ndata_dict, verbose=False):
        self.model = model
        self.verbose = verbose
        self.ndata_dict = ndata_dict

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
        # assume they come in the correct order
        vl_dict = {}
        total_points = 0
        for loss, (exp_name, npoints) in zip(loss_list[1:], self.ndata_dict.items()):
            vl_dict[exp_name] = loss / npoints
            total_points += npoints

        total_loss = loss_list[0] / total_points

        return total_loss, vl_dict

    @property
    def weights(self):
        return self.__weights

    @weights.setter
    def weights(self, weights):
        self.model.set_weights(weights)

    @weights.getter
    def weights(self):
        return self.model.get_weights()

    @property
    def loss(self):
        return self.__loss
    @loss.getter
    def loss(self):
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

    def __init__(self, threshold):
        self.threshold = threshold
        self.positivity = False

    def check_positivity(self, history_object, pos_key="POS"):
        """
            This function receives a history object and look for entries
            with the keyname: POS_something

            A `history_object` is a dictionary of entries in the form
                {'name' : loss}
        """
        positivity_loss = 0.0
        for key, item in history_object.items():
            if key.startswith(pos_key):
                positivity_loss += item[-1]
        if positivity_loss > self.threshold:
            self.positivity = False
        else:
            self.positivity = True
        return self.positivity

    def __call__(self):
        if self.positivity:
            return POS_OK
        else:
            return POS_BAD




