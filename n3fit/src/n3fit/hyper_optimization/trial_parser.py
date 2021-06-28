"""
    Read a folder with trials from a hyperopt run
    and parse it into a number of Trial objects.

    The current algorithm gets the x best trials
    ordered by reward (1/loss)
    and then randomly samples from the array with a weight according
    to the reward.

    Example
    -------
    >>> from n3fit.hyper_optimization.trial_parser import produce_trials
    >>> trials = produce_trials("trial_folder" , size=10)
"""
import json
import logging
from glob import glob
import numpy as np


log = logging.getLogger(__name__)


class Trial:
    """Trial super class"""

    def __init__(self, trial_dict, minimum_losses=4, original_parameters=None):
        self._trial_dict = trial_dict
        self._minimum_losses = minimum_losses
        self._original_params = original_parameters
        self._reward = None

    @property
    def reward(self):
        """Return and cache the reward value"""
        if self._reward is None:
            self._reward = self._evaluate()
        return self._reward

    @property
    def loss(self):
        """Return the loss of the hyperopt dict"""
        return self._trial_dict["result"]["loss"]

    @property
    def params(self):
        """Parameters for the fit"""
        trial_results = self._trial_dict["misc"]["space_vals"]
        if self._original_params:
            hyperparameters = dict(self._original_params, **trial_results)
        else:
            hyperparameters = trial_results
        # Ensure that no hyperparameters has the wrong shape/form
        # 1. For n3fit, activation_per_layer must be a list
        apl = hyperparameters["activation_per_layer"]
        if isinstance(apl, str):
            apl = [apl]*(len(hyperparameters["nodes_per_layer"])-1) + ["linear"]
            hyperparameters["activation_per_layer"] = apl
        hyperparameters["epochs"] = 250
        return hyperparameters

    # Slightly fake a dictionary behaviour with the params attribute
    def __getitem__(self, item):
        return self.params[item]

    def get(self, item, default=None):
        return self.params.get(item, default)
    #####

    def __str__(self):
        strs = ["Parameters:"] + [f"  {i}: {k}" for i, k in self.params.items()]
        strs.append(f"Reward: {self.reward}")
        str_out = "\n".join(strs)
        return str_out

    def _evaluate(self):
        """Evaluate the reward function for a given trial
        Reward = 1.0/loss
        """
        result = self._trial_dict["result"]
        if result.get("status") != "ok":
            return False
        # If we don't have enough validation losses, fail
        val_loss = result["kfold_meta"].get("validation_losses", [])
        if self._minimum_losses and len(val_loss) < self._minimum_losses:
            return False
        return 1.0 / result["loss"]

    def __gt__(self, another_trial):
        """Return true if the current trial is preferred
        when compared to the target"""
        if not another_trial.reward:
            return True
        if not self._reward:
            return False
        return self._reward > another_trial._reward

    def __lt__(self, another_trial):
        """Return true if the target trial is preferred
        when compared to the current (self) one"""
        return not self > another_trial


def get_trials_from_fit(fitfolder, original_parameters=None):
    """Read all trials from a given fit"""
    all_trial_files = glob(f"{fitfolder}/nnfit/*/tries.json")

    all_trials = []
    for trial_file in all_trial_files:
        with open(trial_file, "r") as tf:
            all_trials += [Trial(i, original_parameters=original_parameters) for i in json.load(tf)]
    return all_trials


def produce_trials(folder, original_parameters=None, size=1, sigma=4.0, quick_trim=50):
    """Sample a trial from a hyperopt run
    according to a distribution given by their losses.
    In the limits:
    a sigma of 0 will effectively sample from a random distribution
    and a sigma of infinity will just return the best
    You can limit yourself to the best `quick_trim`.

    By default return just 1 trial (`size=1`) but can ask
    for as many as you want (without replacement)
    """
    all_trials_raw = get_trials_from_fit(folder, original_parameters=original_parameters)
    # Drop the failing trials
    all_trials = list(filter(lambda i: i.reward, all_trials_raw))
    # Get the n with the biggest reward
    best_x = sorted(all_trials)[-quick_trim:]
    if len(best_x) < quick_trim:
        log.warning("Not enough trials to reach %d best, only %d found", quick_trim, len(best_x))
    if len(best_x) < size:
        log.warning("Asked for %d trials, more than what we can get with these options", size)
    # Compute weights proportionally to the reward (which goes from 0 (worst) to 1 (best, loss=1))
    rewards = np.array([i.reward for i in best_x])
    weight_raw = np.exp(sigma * rewards ** 2)
    total = np.sum(weight_raw)
    weights = weight_raw / total
    return np.random.choice(best_x, replace=False, size=size, p=weights)
