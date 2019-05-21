#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Custom hyperopt trial object for persistent file storage
"""
import json
from hyperopt import Trials, space_eval


def space_eval_trial(space, trial):
    for_eval = {}
    for k, v in trial["misc"]["vals"].items():
        if len(v) == 0:
            for_eval[k] = None
        else:
            for_eval[k] = v[0]
    return space_eval(space, for_eval)


class FileTrials(Trials):
    """
    Stores trial results on the fly.
    """

    def __init__(self, replica_path, log=None, parameters=None, exp_key=None, refresh=True):
        self._store_trial = False
        self._json_file = "{0}/tries.json".format(replica_path)
        self._parameters = parameters
        if log:
            self.log_info = log.info
        else:
            self.log_info = print
        super(FileTrials, self).__init__(exp_key=exp_key, refresh=refresh)

    def refresh(self):
        super(FileTrials, self).refresh()

        # write json to disk
        if self._store_trial:
            self.log_info(f"Storing scan in {self._json_file}")
            local_trials = []
            for idx, t in enumerate(self._dynamic_trials):
                local_trials.append(t)
                local_trials[idx]["misc"]["space_vals"] = space_eval_trial(self._parameters, t)

            all_to_str = json.dumps(local_trials, default=str)
            with open(self._json_file, "w") as f:
                f.write(all_to_str)

    def new_trial_ids(self, N):
        self._store_trial = False
        return super(FileTrials, self).new_trial_ids(N)

    def new_trial_docs(self, tids, specs, results, miscs):
        self._store_trial = True
        return super(FileTrials, self).new_trial_docs(tids, specs, results, miscs)
