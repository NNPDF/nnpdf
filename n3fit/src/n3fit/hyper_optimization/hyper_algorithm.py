"""
    This module contains functions dedicated to process the json dictionaries
"""

import itertools
import numpy as np

fail_threshold = 75
fail_reward = -100
key_good = "good"
key_loss = "loss"

def compute_reward(mdict, biggest_ntotal):
    """
    Given a combination dictionary computes the reward function:

    If the fail rate for this combination is above the fail threshold, rewards is -100

    The formula for the reward function is:
        rwd = combination_weight*( (0.5-fail_rate) - dispersion/good_rate )
    """
    # Check the fail rate
    fail_rate = mdict["fail_rate"]
    if fail_rate > fail_threshold:
        return fail_reward
    # Compute the good rate
    n_good = mdict["n_good"]
    n_total = mdict["n_total"]
    rate_good = n_good / n_total * 100.0
    # Compute the dispersion
    dispersion = mdict["median"] - mdict["best_loss"]
    # And with all of this, compute a reward function
    reward = (50.0 - fail_rate)/100.0
    reward -= dispersion/rate_good # The less spread the better
    # Finally scale it with the total number of points
    combination_weight = 1.0 + n_total/biggest_ntotal
    return reward*combination_weight

def parse_keys(dataframe, keys):
    """
    Receives a dataframe and a set of keys
    Looks into the dataframe to read the possible values of the keys

    Returns a dictionary { 'key' : [possible values] },

    # Arguments:
        - `dataframe`: a pandas dataframe
        - `keys`: keys to combine

    # Returns:
        - `key_info`: a dictionary with the possible values for each key
    """
    key_info = {}
    for key_name in keys:
        # Get all values
        all_possible_values = dataframe[key_name]
        # Remove nans and duplicates
        removed_duplicates = sorted(list(set(all_possible_values)))
        removed_nans = filter(lambda x: x == x, removed_duplicates)
        # If there's anything left, add it to the dictionary
        if removed_nans:
            key_info[key_name] = removed_nans
    return key_info

def get_combinations(key_info, n):
    """
    Get all combinations of n elements for the keys and values given in the dictionary key_info:
    { 'key' : [possible values] }, if key_info is None, the dictionary used to initialize the class will be used
    Returns a list of dictionaries such that:
    [
    {'key1' : val1-1, 'key2', val2-1 ... },
    {'key1' : val1-2, 'key2', val2-1 ... },
    ...
    ]

    # Arguments:
        - `key_info`: dictionary with the possible values for each key
        - `n`: elements to combine

    # Returns:
        - `all_combinations`: A list of dictionaries of parameters
    """
    # If we don't have enough keys to produce n combinations, return empty
    if len(key_info) < n:
        return []
    # First generate the combinations of keys
    key_combinations = itertools.combinations(key_info, n)
    all_combinations = []
    # Now, for each combination of keys we have to give values
    for keys in key_combinations:
        # Generate a list of tuples with the values of the keys
        # i.e., something like [ (val1-1, val1-2, val1-3...), (val2-1, val2-2...) ... ]
        list_of_items = [key_info[key] for key in keys]
        # Now combine all these items, which is what we actually want
        items_combinations = itertools.product(*list_of_items)
        # Now we want to put things back in the form of a dictionary of parameters
        for values in items_combinations:
            # Each values comes in the same order as `keys`
            new_dictionary = dict(zip(keys, values))
            all_combinations.append(new_dictionary)
    return all_combinations

def get_slice(dataframe, query_dict):
    """
    Returns a slice of the dataframe where some keys match some values
    keys_info must be a dictionary {key1 : value1, key2, value2 ...}
    # Arguments:
        - `dataframe`: a pandas dataframe
        - `query_dict`: a dictionary of combination as given by `get_combinations`
    """
    df_slice = dataframe
    for key, value in query_dict.items():
        key_column = df_slice[key]
        # Check whether for this slice the values are not NaN
        if len(key_column.dropna()) == 0 and len(key_column) != 0:
            return None
        # TODO with continuous keys we don't have specific values!
        df_slice = df_slice[key_column == value]
    return df_slice

def process_slice(df_slice):
    """
    Function to process a slice into a dictionary with useful stats
    If the slice is None it means the combination does not apply

    # Arguments:
        - `df_slice`: a slice of a pandas dataframe

    # Returns:
        - `proc_dict`: a dictionary of stats
    """
    # First check whether there's anything inside the slice
    n_total = len(df_slice)
    if df_slice is None or n_total == 0:
        proc_dict = {"skip" : True}
        return proc_dict
    else:
        proc_dict = {"skip" : False}
    # Get the good values
    good = df_slice[df_slice[key_good]]
    # Get raw stats
    n_good = len(good)
    n_failed = n_total - n_good
    fail_rate = n_failed / n_total * 100.0
    # Now get the distribution of the (good) losses
    good_losses = good[key_loss]
    best_loss = good_losses.min()
    std_dev = good_losses.std()
    median = good_losses.median()
    # Fill the dictionary
    proc_dict["n_failed"] = n_failed
    proc_dict["n_good"] = n_good
    proc_dict["n_total"] = n_total
    proc_dict["fail_rate"] = fail_rate
    proc_dict["best_loss"] = best_loss
    proc_dict["std"] = std_dev
    proc_dict["median"] = median
    return proc_dict


def study_combination(dataframe, query_dict):
    """
    Given a dataframe and a dictionary of {key1 : value1, key2: value2}
    returns a dictionary with a number of stats for that combination

    # Arguments:
        - `dataframe`: a pandas dataframe
        - `query_dict`: a dictionary for a combination as given by `get_combinations`

    # Returns:
        - `proc_dict`: a dictionary of the "statistics" for this combination
    """
    # Get the slice corresponding to this combination
    df_slice = get_slice(dataframe, query_dict)
    proc_dict = process_slice(df_slice)
    proc_dict["slice"] = df_slice
    proc_dict["combination"] = query_dict
    return proc_dict

def dataframe_removal(dataframe, hit_list):
    """
    Removes all combinations defined in hit_list from the dataframe.
    The hit list is list of dictionaries containing the 'slice' key
    where 'slice' must be a slice of 'dataframe'

    # Arguments:
        - `dataframe`: a pandas dataframe
        - `hit_list`: the list of element to remove

    # Returns:
        - `new_dataframe`: the same dataframe with all elements from hit_list removed
    """
    if not hit_list:
        return dataframe
    # I think I am failing to understand how the index object works
    # this looks unnecesaryly verbose
    # I am just getting all the indices from all the different combinations
    # making sure there are no duplicates and then dropping them from the dataframe
    indices_to_remove = hit_list[0]["slice"].index
    for victim in hit_list[1:]:
        indices_to_remove = indices_to_remove.append(victim["slice"].index)
    indices_to_drop = indices_to_remove.drop_duplicates()
    new_dataframe = dataframe.drop(indices_to_drop)
    return new_dataframe


def autofilter_dataframe(dataframe, keys,  n_to_combine = 2, n_to_kill = 2, threshold=-50):
    """
    Receives a dataframe and a list of keys.
    Creates combinations of `n_to_combine` keys and computes the reward
    Finally removes from the dataframe the `n_to_kill` worse combinations

    Anything under `threshold` will be removed and will not count towards the `n_to_kill`
    (by default `threshold` = -50 so only things which are really bad will be removed)

    # Arguments:
        - `dataframe`: a pandas dataframe
        - `keys`: keys to combine
        - `n_to_combine`: how many keys do we want to combine
        - `n_to_kill`: how many combinations to kill
        - `threshold`: anything under this reward will be removed

    # Returns:
        - `dataframe_sliced`: a slice of the dataframe with the weakest combinations
                              removed
    """
    # Step 0: read from the dataframe the content of the keys
    key_info = parse_keys(dataframe, keys)
    # Step 1: get the combinations
    combinations = get_combinations(key_info, n_to_combine)
    # Step 2: run through all possible combinations and compute stats
    result_list = []
    hit_list = []
    biggest_ntotal = 1
    for combination in combinations:
        processed_dict = study_combination(dataframe, combination)
        if processed_dict["skip"]:
            continue
        if processed_dict["n_total"] > biggest_ntotal:
            biggest_ntotal = processed_dict["n_total"]
        result_list.append(processed_dict)
    # Step 3: compute reward
    for processed_dict in result_list:
        reward = compute_reward(processed_dict, biggest_ntotal)
        if reward < fail_reward:
            hit_list.append(combination)
        processed_dict["reward"] = reward
    # Step 4: Order the results by reward
    result_list.sort(key = lambda i: i["reward"])
    # Step 5: Add the n-last to the list of combinations to remove
    n_removal = len(hit_list) + n_to_kill
    hit_list += result_list[-n_removal:]
    for i in hit_list:
        print("Removing {0} with reward {1}".format(i['combination'], i['reward']))
    # Step 6: remove the bad guys from the dataframe
    new_dataframe = dataframe_removal(dataframe, hit_list)
    return new_dataframe





# code graveyard, to be removed in a few commits

import itertools
import numpy as np


class HyperAlgorithm:
    def __init__(
        self,
        key_info,
        verbose=True,
        k_good="good",
        k_loss="loss",
        threshold_reward=0.0,
        fail_threshold=75,
        how_many_to_kill=2,
        remove_by_key=False,
    ):
        """
            key_info: dictionary of {'keyname' : [list,of,possible,values]}
        """

        self.key_info = key_info

        self.k_good = k_good
        self.k_loss = _loss

        self.verbose = verbose
        self.remove_by_key = remove_by_key

        self.threshold_reward = threshold_reward
        self.failing_reward = -100
        self.fail_threshold = fail_threshold

        self.how_many_to_kill = how_many_to_kill

        # Some variables are not discrete, find them
        update_dict = {}
        for key, values in key_info.items():
            try:  # Remove nan if there are
                np_values = np.array(values)
                np_values = np_values[~np.isnan(np_values)]
            except:
                np_values = values
            l = len(np_values)
            if l > 10:  # break it into five chuncks
                np_values.sort()
                new_vals = []
                step = int(l / 5)
                for i in range(0, l - step, step):
                    ini = np_values[i]
                    fin = np_values[i + step]
                    new_vals.append((ini, fin))
                update_dict[key] = new_vals
        key_info.update(update_dict)

        self.continuous = update_dict
        self.biggest_total = 1
        self.perfect_reward = 100.0

    def __dataframe_killer(self, dataframe, hit_list):
        """
        Remove from the dataframe all entries included in the hit_list
        """
        if not hit_list:
            return dataframe
        indices_to_remove = hit_list[0]["slice"].index
        for victim in hit_list[1:]:
            indices_to_remove = indices_to_remove.append(victim["slice"].index)
        indices_to_drop = indices_to_remove.drop_duplicates()
        return dataframe.drop(indices_to_drop)

    def __verbose_elimination(self, md, printout=True):
        """
        Show a log of everything that is being eliminated
        if there is only one item in the dictionary,
        it also removes it from the combination list just to avoid looking at useless combinations
        """
        combination = md["combination"]
        str_list = []
        for key, item in combination.items():
            str_list.append("{0}={1}".format(key, item))
        reward = md["reward"]
        if printout:
            self.__print("Removing: {0}\n  reward: {1:.5f}".format("\n          ".join(str_list), reward))

    def __print(self, string, **kwargs):
        if self.verbose:
            print(string, **kwargs)

    def __compute_reward(self, md):
        """ Computes a reward function from the information we have """
        frate = md["f_rate"]
        ngood = md["n_good"]
        ntotal = md["n_total"]
        # If the fail rate is greater than 75%, kill it
        if frate > self.fail_threshold:
            return self.failing_reward
        dstd = md["median"] - md["best_loss"]
        rate_good = ngood / ntotal * 100.0

        # The reward can be negative
        reward = (50.0 - frate) / 100.0
        reward -= dstd / rate_good  # the less spread the better
        # Finally scale it with the total number of points
        weight_points = 1.0 + ntotal / self.biggest_total
        reward *= weight_points

        return reward

    def print_ranges(self, dataframe, keys=None, dataframe_original=None):
        """
        Receives a dataframe and a number of keys and prints the surviving keys and ranges
        If dataframe_original is given, will also print (for discrete values)
        how many of each have survived from the total
        """
        # First get the keys we are printing
        keys = self.key_info.keys()

        # Now compute all rewards
        key_rewards = {}
        for key in keys:
            combinations = self.get_combinations(1, single_key=key)
            key_rewards[key] = []
            for combination in combinations:
                md = self.study_combination(dataframe, combination)
                key_rewards[key].append(md)

        if dataframe_original is None:
            prev = ""
        else:
            prev = "of {0}".format(len(dataframe_original))
        print("The total number of survivors is {0} {1}".format(len(dataframe), prev))

        # Now do the printing
        for key, results in key_rewards.items():
            print("For key: {0}".format(key))
            list_str = []
            results.sort(key=lambda i: -i["reward"])
            for md in results:
                val = md["combination"][key]
                reward = md["reward"]
                list_str.append("         Reward: {1:.3f}: {0} ".format(val, reward))
            print("\n".join(list_str))
        return

        for key in keys:
            df_slice = dataframe[key]
            values = set(df_slice)
            if len(values) > 5:
                ini = df_slice.min()
                fin = df_slice.max()
                string = " range: ({0}, {1}) ".format(ini, fin)
            else:
                if dataframe_original is not None:
                    originals = dataframe_original[key].value_counts()
                    news = df_slice.value_counts()
                    string = " values: "
                    for val in values:
                        string += "{0} ({1}/{2} {3:2.0f}%), ".format(
                            val, news[val], originals[val], news[val] / originals[val] * 100.0
                        )
                    string = string[:-2]
                else:
                    string = " values: {0}".format(values)

            self.__print("Key: {0}\n > {1}".format(key, string))

    def __process_slice(self, df_slice):
        """ Function to process a slice into a dictionary with interesting stats
        If the slice is None it means the combination does not apply
        """
        if df_slice is None:
            md = {"skip": True}
            return md
        else:
            md = {"skip": False}

        good = df_slice[df_slice[self.k_good]]
        md["n_total"] = len(df_slice)
        md["n_good"] = len(good)
        if md["n_good"] > self.biggest_total:
            self.biggest_total = md["n_good"]
        md["n_failed"] = md["n_total"] - md["n_good"]
        if md["n_total"] == 0:
            md["skip"] = True
            return md
        else:
            f_rate = md["n_failed"] / md["n_total"] * 100
        md["f_rate"] = f_rate
        # Compute for each the best loss (among the good ones) and the std dev
        good_losses = good[self.k_loss]
        md["best_loss"] = good_losses.min()
        md["std"] = good_losses.std()
        md["median"] = good_losses.median()
        md["reward"] = self.__compute_reward(md)
        return md

    def get_slice(self, dataframe, keys_info):
        """
        Returns a slice of the dataframe where some keys match some values
        keys_info must be a dictionary {key1 : value1, key2, value2 ...}
        """
        df_slice = dataframe
        for key, value in keys_info.items():
            # First check whether for this slice the values for this key are something different from NaN
            key_column = df_slice[key]
            if len(key_column.dropna()) == 0 and len(key_column) != 0:
                return None
            if key in self.continuous.keys():
                minim = value[0]
                maxim = value[1]
                df_slice = df_slice[key_column >= minim]
                df_slice = df_slice[key_column <= maxim]
            else:
                df_slice = df_slice[key_column == value]
        return df_slice

    def get_combinations(self, n, key_info=None, single_key=None):
        """
        Get all combinations of n elements for the keys and values given in the dictionary key_info:
        { 'key' : [possible values] }, if key_info is None, the dictionary used to initialize the class will be used
        Returns a list of dictionaries such that:
        [
        {'key1' : val1-1, 'key2', val2-1 ... },
        {'key1' : val1-2, 'key2', val2-1 ... },
        ...
        ]
        """
        if key_info is None:
            key_info = self.key_info
        if len(key_info) < n:
            return []

        if single_key:
            single_options = key_info[single_key]
            key_info = {single_key: single_options}

        key_combinations = itertools.combinations(key_info, n)
        all_combinations = []
        for keys in key_combinations:
            list_of_items = [key_info[key] for key in keys]
            # list_of_items is a list of possible values: [ (val1-1, val1-2, val1-3...), (val2-1, val2-2...) ... ]
            items_combinations = itertools.product(*list_of_items)
            for values in items_combinations:
                nd = dict(zip(keys, values))
                all_combinations.append(nd)

        return all_combinations

    def study_combination(self, dataframe, combination):
        """
        Given a dataframe and a dictionary of {key1 : value1, key2: value2}
        returns a dictionary with a number of stats for that combination
        """
        df_slice = self.get_slice(dataframe, combination)
        # If the combination does not applies (i.e., NaN) just skip it
        md = self.__process_slice(df_slice)
        md["slice"] = df_slice
        md["combination"] = combination
        return md

    def __sort_to_survive(self, info_dict, how_many_to_save):
        """
        Receives a list of dictionaries with an entry called reward, returns the
        'n' best
        """
        sorted_by_reward = sorted(info_dict, key=lambda i: i["reward"])
        sorted_by_reward.reverse()
        # Everyone with perfect reward is saved
        save_list = []
        for i in sorted_by_reward:
            if how_many_to_save == 0:
                break
            if i["reward"] == self.perfect_reward:
                save_list.append(i)
            else:
                save_list.append(i)
                how_many_to_save -= 1
        return save_list

    def __sort_to_kill(self, info_dict, by_key=False, how_many_to_kill=None):
        """
        Receives a list of dictionaries with an entry called reward, returns the
        'n' worst ones.
        If by key is True, returns the 'n' worst for each set of keys
        """
        if how_many_to_kill is None:
            how_many_to_kill = self.how_many_to_kill
        if by_key:
            list_by_key = {}
            for md in info_dict:
                keys = md["keys"]
                if keys in list_by_key.keys():
                    list_by_key[keys].append(md)
                else:
                    list_by_key[keys] = [md]
            hit_list = []
            for _, items in list_by_key.items():
                hit_list += self.__sort_to_kill(items)
        else:
            sorted_by_reward = sorted(info_dict, key=lambda i: i["reward"])
            hit_list = sorted_by_reward[:how_many_to_kill]
        return hit_list

    def remove_failing(self, dataframe, n_combinations=1, key=None, fail_threshold=None, verbose=True):
        """
        Remove entries with failing reward from the dataframe
        """
        if fail_threshold is None:
            set_fail_threshold = self.failing_reward
        elif fail_threshold == "median":
            set_fail_threshold = 100
        else:
            set_fail_threshold = fail_threshold
        combinations = self.get_combinations(n_combinations, single_key=key)
        hit_list = []
        all_combinations = []

        for combination in combinations:
            md = self.study_combination(dataframe, combination)
            if md["skip"]:
                continue
            reward = md["reward"]
            if fail_threshold == "median":
                all_combinations.append(md)
            elif reward <= set_fail_threshold:
                hit_list.append(md)

        if fail_threshold == "median":
            sorted_by_reward = sorted(all_combinations, key=lambda i: i["reward"])
            # do like thanos
            half = int(len(sorted_by_reward) / 2)
            hit_list = sorted_by_reward[:half]

        new_dataframe = self.__dataframe_killer(dataframe, hit_list)
        if verbose:
            for i in hit_list:
                print("Eliminated {0} with reward: {1}".format(i["combination"], i["reward"]))
        return new_dataframe

    def check_wrapper_weaklings_killer(self, dataframe, n_combinations=1, key=None, verbose=False, kill_threshold=None):
        """
        It starts killing the entry with the worst reward and won't stop until all remaining entries
        for the given number of combinations, have a positive reward
        """
        if kill_threshold is None:
            set_kill_threshold = self.threshold_reward
        elif kill_threshold is "median":
            set_kill_threshold = 100  # ensures that someone dies here
        else:
            set_kill_threshold = kill_threshold

        # 1 - Get combinations
        combinations = self.get_combinations(n_combinations, single_key=key)
        # 2 - Run through all possible combinations and kill the weakest one
        worst_reward = 1e5
        victim = None
        all_combinations = []
        for combination in combinations:
            md = self.study_combination(dataframe, combination)
            if md["skip"]:
                continue
            all_combinations.append(md)
            reward = md["reward"]
            if reward < worst_reward and reward <= set_kill_threshold:
                victim = md
                worst_reward = reward
        if verbose:
            all_combinations.sort(key=lambda i: i["reward"])
            for i in all_combinations:
                print("{0}, {1}".format(i["reward"], i["combination"]))

        if kill_threshold == "median":  # Set the median which will be necessary for the rest of the iterations
            rewards = [i["reward"] for i in all_combinations]
            kill_threshold = np.median(rewards)
            print(" > > Killing everything with reward below {0}".format(kill_threshold))

        if victim is None:
            # If there are no points below the reward threshold, go out and return the remaining dataframe
            return dataframe

        # 3 - Kill the victim(s)
        new_dataframe = self.__dataframe_killer(dataframe, [victim])
        self.__verbose_elimination(victim)

        # 4 - Keep on murdering
        return self.check_wrapper_weaklings_killer(
            new_dataframe, n_combinations, key=key, kill_threshold=kill_threshold
        )

    def check_wrapper(self, dataframe, n_combinations, kill=False, save_n_best=None):
        """
        Input a dataframe and the number of possible combinations and prints out nice statistics
        """

        def get_str(comb, md):
            key_strs = ["{0}={1}".format(key, item) for key, item in comb.items()]
            comb_str = ", ".join(key_strs)
            f_rate = md["f_rate"]
            reward = md["reward"]
            n_trials = md["n_total"]
            failure_str = "Fail: {0:2.3f}% Reward: {2:0.3f} ({1:3d} trials)".format(f_rate, n_trials, reward)
            return comb_str, failure_str, f_rate

        self.__print("\n > Check combinations of {0} keys".format(n_combinations))

        # 1 - Get all possible combinations of the keys
        combinations = self.get_combinations(n_combinations)
        lines_to_print = []
        siz = 0
        # 2 - Now run for all possible combinations and parse the information
        hit_list = []
        info_dict = []
        for combination in combinations:
            md = self.study_combination(dataframe, combination)
            if md["skip"]:
                continue
            reward = md["reward"]
            if reward <= self.threshold_reward:
                # Don't bother with this point, kill it later
                hit_list.append(md)
            else:
                md["keys"] = tuple(sorted((combination.keys())))
                info_dict.append(md)
                comb_str, fail_str, _ = get_str(combination, md)
                if len(comb_str) > siz:
                    siz = len(comb_str)
                lines_to_print.append((comb_str, fail_str, reward))

        # 3 - If kill == true, remove the slices in the hit list
        if kill:
            # 3.1 - Sort the items by reward and get the 'how_many_to_kill' with the smallest reward
            hit_list += self.__sort_to_kill(info_dict, by_key=self.remove_by_key)
            print("Removing the worst combinations:")
            new_dataframe = self.__dataframe_killer(dataframe, hit_list)
        else:
            new_dataframe = dataframe

        # 4 - If save_n_best, do the opposite
        if save_n_best:
            survivors_list = self.__sort_to_survive(info_dict, save_n_best)
            # Everyone with a 100.0 reward must survive
            indices_to_get = None
            for survivor in survivors_list:
                if indices_to_get is not None:
                    indices_to_get = indices_to_get.append(survivor["slice"].index)
                else:
                    indices_to_get = survivor["slice"].index
            indices_to_get = indices_to_get.drop_duplicates()
            new_dataframe = dataframe.ix[indices_to_get]

        # 5 - Now print all lines, making sure the align nicely
        print("\n")
        lines_to_print.sort(key=lambda x: x[-1])
        # If saving only the best X, print only those
        if save_n_best:
            lines_to_print = lines_to_print[-save_n_best:]
        base_line = "{0:" + str(siz) + "s} {1}"
        for line in lines_to_print:
            self.__print(base_line.format(line[0], line[1]))

        return new_dataframe
