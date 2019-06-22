import numpy as np
from matplotlib import pyplot as pl


class Stat_Info:
    def __init__(
        self,
        vl_model,
        vl_data,
        ndata_dict,
        threshold_positivity=1e-6,
        chi2_threshold=3.0,
        save_all_each=100,
        total_epochs=0,
        stopping_epochs=7000,
        dont_stop=False,
    ):
        """
        To instantiate the Stats_Info object we need:
            - vl_model: The validation  model, which has the same PDF layer as the training model
            - vl_data: The validation exp data to compare to
            - ndata_dict: A dictionary with the number of points of each experiment in the format:
                exp: npoints in the training model for experiment exp
                exp_vl : npoints in the training model for experiment exp
            - chi2_threshold: at which training chi do we start checking for a stopping point?
            - save_all_each: in which interval should we save all chi2 per experiment?
            - stopping_epochs: after a growth in chi2 occurs, we activate the stopping monitor
                if the chi2 doesnt go below the minimum after stopping_epochs, we stop
            - dont_stop: dont_stop anyway

        """  # TODO: this function has outgrown its initial scope and some refactoring is needed

        self.vl_model = vl_model
        self.vl_data = vl_data

        self.stop_now = False
        self.nan_found = False
        self.stopping_epochs = stopping_epochs
        self.stopping_degree = 0
        self.total_epochs = total_epochs
        self.epoch_of_the_stop = total_epochs
        self.dont_stop = dont_stop
        self.total_epochs = total_epochs
        self.not_reloaded = True
        self.save_all_each = save_all_each

        # Threshold for positivity
        self.threshold_positivity = threshold_positivity

        # Best model stats
        self.threshold = chi2_threshold

        # Number of points per experiment
        self.ndata_tr = max(ndata_dict["total_tr"], 1)
        self.ndata_vl = max(ndata_dict["total_vl"], 1)
        self.ndata_dict = ndata_dict
        self.exp_names = []
        for i in ndata_dict.keys():
            if i.endswith("_vl") or i.startswith(("POS", "total", "ALL")):
                continue
            self.exp_names.append(i)

        self._clean_best_model()
        self._clean_full_stats()

    def _clean_best_model(self):
        self.best_chi2 = 1e9
        self.e_best_chi2 = 0
        self.w_best_chi2 = None
        self.positivity = False
        self.saved_history = None

    def _clean_full_stats(self):
        # Complete stats
        self.epochs = []
        self.tr_chi2 = []
        self.vl_chi2 = []
        self.all_tr_chi2 = {}
        self.all_vl_chi2 = {}
        self.all_epochs = []
        for exp in self.exp_names:
            self.all_tr_chi2[exp] = []
            self.all_vl_chi2[exp] = []
        self.all_tr_chi2["Total"] = []
        self.all_vl_chi2["Total"] = []

    def reset(self, stopping_epochs=None):
        """
        Sometimes it might be necessary to reset the statistics
        This only require resetting some of the variables and not all
        """
        self._clean_full_stats()
        self._clean_best_model()
        self.stop_now = False
        self.epoch_of_the_stop = self.total_epochs
        self.not_reloaded = True
        self.stopping_degree = 0
        if stopping_epochs:
            self.stopping_epochs = stopping_epochs

    def _compute_validation(self):
        """ Evaluates the validation model
        vl_list: a list of chi2, where the first element is the sum of all others
        If there is only one experiment/dataset, vl_list is a float instead

        Returns the list already normalised to the number of points
        """
        vl_list = self.vl_model.evaluate(x=None, y=self.vl_data, verbose=False, batch_size=1)
        vl_dict = {}
        try:
            total = vl_list[0] / self.ndata_vl
            vl_list = vl_list[1:]
        except IndexError:
            total = vl_list / self.ndata_vl
            vl_list = []

        # If validation dont include all experiments, there is something wrong in the modle
        # and since evaluate returns a list with no names,
        # it is impossible to know which experiment corresponds to which chi2
        if len(vl_list) != len(self.exp_names):
            vl_list = [0.0 for i in self.exp_names]

        # This loop relies on the information that comes through the input dictionaries to be accurate
        # It also relies in the ordered dictionaries introduced at some point in python 3.6/7
        for chi2, exp in zip(vl_list, self.exp_names):
            vl_dict[exp] = chi2 / self.ndata_dict[exp]

        return total, vl_dict

    def _parse_training(self, training_info):
        """ Receives an object containing the training chi2 and parses it
        in the form of a { 'experiment' : chi2 } dictionary, where chi2
        is already normalised to the number of points per experiment """
        try:
            hobj = training_info.history
        except AttributeError:  # So it works whether we pass the out our the out.history
            hobj = training_info

        # Since the history object can have more than one loss (if epochs != 1) but it is always a list
        # Save them as the mean

        total = np.mean(hobj["loss"]) / self.ndata_tr
        tr_chi2 = {}
        for exp in self.exp_names:
            chi2 = np.mean(hobj[exp + "_loss"]) / self.ndata_dict[exp]
            tr_chi2[exp] = chi2

        return total, tr_chi2

    def monitor_chi2(self, training_info, epoch):
        """ Function to be called at the end of every epoch.
        Stores the total chi2 of the training set as well as the
        total chi2 of the validation set.
        If the training chi2 is below a certain threshold,
        stores the state of the model which gave the minimum chi2
        as well as the epoch in which occurred
        If the epoch is a multiple of save_all_each then we also save the per-exp chi2

        Returns True if the run seems ok and False if a NaN is found
        """
        vl_chi2, all_vl = self._compute_validation()
        tr_chi2, all_tr = self._parse_training(training_info)

        if np.isnan(tr_chi2):
            print(" > NaN found, stopping activated")
            self.stop_now = True
            self.nan_found = True
            self.epoch_of_the_stop = epoch
            self.reload_model()
            if not self.w_best_chi2:
                self.best_chi2 = 1e10
            return False

        if vl_chi2 == 0.0:  # If there is no masking, training and validation are the same thing
            vl_chi2 = tr_chi2
            all_vl = all_tr

        self.epochs.append(epoch)
        self.tr_chi2.append(tr_chi2)
        self.vl_chi2.append(vl_chi2)

        if tr_chi2 < self.threshold:
            self.stopping_degree += 1

            if vl_chi2 < self.best_chi2 and self.check_positivity(training_info):
                self.positivity = True
                self.best_chi2 = vl_chi2
                self.e_best_chi2 = epoch
                self.w_best_chi2 = self.vl_model.get_weights()
                self.stopping_degree = 0
                self.saved_history = training_info
            elif not self.w_best_chi2:
                self.stopping_degree = 0

        if (epoch + 1) % 100 == 0:
            self.all_epochs.append(epoch + 1)
            for exp in self.exp_names:
                self.all_tr_chi2[exp].append(all_tr[exp])
                self.all_vl_chi2[exp].append(all_vl[exp])
            self.all_tr_chi2["Total"].append(tr_chi2)
            self.all_vl_chi2["Total"].append(vl_chi2)
            self.print_current_stats()

        if self.stopping_degree > self.stopping_epochs:
            self.stop_now = True
            self.epoch_of_the_stop = epoch
            self.reload_model(force=True)

        return True

    def reload_model(self, force=False):
        if (self.not_reloaded or force) and self.w_best_chi2:
            self.not_reloaded = False
            self.vl_model.set_weights(self.w_best_chi2)

    def good_stop(self):
        """ Returns true if and only if the stopping criteria was actually fulfilled
        and we didnt stop by error/nan """
        return self.positivity and not self.nan_found

    def stop_here(self):
        if self.dont_stop:
            return False
        else:
            return self.stop_now

    def loss(self):
        return self.best_chi2

    def tr_loss(self):
        i = self.epochs.index(self.e_best_chi2)
        return self.tr_chi2[i]

    def print_current_stats(self):
        """
        Prints the last validation and training loss saved
        """
        epoch = self.all_epochs[-1]
        epochs = self.total_epochs
        total = self.all_tr_chi2["Total"][-1]
        total_str = "\nAt epoch {0}/{1}, total loss: {2}".format(epoch, epochs, total)

        partials = ""
        for exp in self.exp_names:
            partials += "{0}: {1:.3f}  ".format(exp, self.all_tr_chi2[exp][-1])

        print(total_str)
        if partials:
            print(" > Partial losses: " + partials)

        total_vl = self.all_vl_chi2["Total"][-1]
        print(" > Validation loss at this point: {0}".format(total_vl))

    def plot(self):
        print("Plotting validation loss")
        pl.plot(self.epochs, self.tr_chi2, label="Training")
        pl.plot(self.epochs, self.vl_chi2, label="Validation")
        pl.legend()
        pl.show()

    def check_positivity(self, training_info):
        """ Checks whether positivity passes """
        chi2 = 0.0
        try:
            hobj = training_info.history
        except AttributeError:  # So it works whether we pass the out our the out.history
            hobj = training_info

        for key, item in hobj.items():
            if "POS" in key:
                if item[-1] > 0.0:
                    chi2 += item[-1]
        if chi2 > self.threshold_positivity:
            return False
        else:
            return True

    # Str functions
    def positivity_pass(self):
        if self.positivity:
            return "POS_PASS"
        else:
            return "POS_VETO"

    def preproc_str(self):
        preprocessing = self.vl_model.get_layer("pdf_prepro")
        alphabeta = preprocessing.get_weights()
        file_list = ["Alpha Beta"]
        for i in range(0, len(alphabeta), 2):
            file_list.append("\n{0} {1}".format(alphabeta[i][0], alphabeta[i + 1][0]))
        file_list.append("\n")
        return file_list

    def chi2exps_str(self):
        file_list = []
        nexp = len(self.exp_names)
        for i, epoch in enumerate(self.all_epochs):
            strout = "\nEpoch: {0} nexp {1}".format(epoch, nexp)
            for exp in self.exp_names:
                tr = self.all_tr_chi2[exp][i]
                vl = self.all_vl_chi2[exp][i]
                strout += "\n{0}: {1} {2}".format(exp, tr, vl)
            tr = self.all_tr_chi2["Total"][i]
            vl = self.all_vl_chi2["Total"][i]
            strout += "\nTotal: training = {0} validation = {1}\n".format(tr, vl)
            file_list.append(strout)
        return file_list
