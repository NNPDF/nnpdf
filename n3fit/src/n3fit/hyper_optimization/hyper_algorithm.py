"""
    This module contains functions dedicated to process the json dictionaries
"""
import itertools
import numpy as np
from collections import namedtuple

MinMax = namedtuple("MinMax", "min max")


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
    reward = (50.0 - fail_rate) / 100.0
    reward -= dispersion / rate_good  # The less spread the better
    # Finally scale it with the total number of points
    combination_weight = 1.0 + n_total / biggest_ntotal
    return reward * combination_weight


def bin_generator(values, max_n=10):
    """
    Receives a list of values. If there are more than `max_n` of them
    and they are numeric, create `max_n` bins.
    If they are already discrete values or there are less than `max_n` options,
    output the same input

    # Arguments:
        - `values`: list of possible values
        - `maximum`: maximum number of allowed different values

    # Returns:
        - `new_vals`: list of tuples with (initial, end) value of the bin
    """
    lval = len(values)
    if lval <= max_n:
        return values
    if not all(isinstance(i, (int, float)) for i in values):
        return values
    values.sort()
    new_vals = []
    step = int(lval / max_n)
    for i in range(0, lval - step, step):
        ini = values[i]
        fin = values[i + step]
        new_vals.append(MinMax(ini, fin))
    return new_vals


def parse_keys(dataframe, keys):
    """
    Receives a dataframe and a set of keys
    Looks into the dataframe to read the possible values of the keys

    Returns a dictionary { 'key' : [possible values] },

    If the values are not discrete then we need to bin it
    let's do this for anything with two many numerical values

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
        removed_duplicates = set(all_possible_values)
        removed_nans = sorted(list(filter(lambda x: x == x, removed_duplicates)))
        # If there's anything left, add it to the dictionary
        if removed_nans:
            # But bin it first in case we find a continous variable
            key_info[key_name] = bin_generator(removed_nans)
    return key_info


def get_combinations(key_info, ncomb):
    """
    Get all combinations of ncomb elements for the keys and values given in the dictionary key_info:
    { 'key' : [possible values] }, if key_info is None, the dictionary used to initialize the class will be used
    Returns a list of dictionaries such that:
    [
    {'key1' : val1-1, 'key2', val2-1 ... },
    {'key1' : val1-2, 'key2', val2-1 ... },
    ...
    ]

    # Arguments:
        - `key_info`: dictionary with the possible values for each key
        - `ncomb`: elements to combine

    # Returns:
        - `all_combinations`: A list of dictionaries of parameters
    """
    # If we don't have enough keys to produce n combinations, return empty
    if len(key_info) < ncomb:
        return []
    # First generate the combinations of keys
    key_combinations = itertools.combinations(key_info, ncomb)
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
        # We need to act differently in the case of continous values we discretized before
        # The way we have to check whether something was continous is to check whether the value
        # is an instance of the class MinMax we created at the top of this file
        if isinstance(value, MinMax):
            minim = value.min
            maxim = value.max
            tmp = df_slice[key_column >= minim]
            df_slice = tmp[key_column <= maxim]
        else:
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
        proc_dict = {"skip": True}
        return proc_dict
    else:
        proc_dict = {"skip": False}
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


def autofilter_dataframe(dataframe, keys, n_to_combine=2, n_to_kill=2, threshold=-50):
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
    biggest_ntotal = 1
    for combination in combinations:
        processed_dict = study_combination(dataframe, combination)
        if processed_dict["skip"]:
            continue
        if processed_dict["n_total"] > biggest_ntotal:
            biggest_ntotal = processed_dict["n_total"]
        result_list.append(processed_dict)
    # Step 3: compute reward
    n_to_remove = n_to_kill
    for processed_dict in result_list:
        reward = compute_reward(processed_dict, biggest_ntotal)
        if reward <= threshold:
            n_to_remove += 1
        processed_dict["reward"] = reward
    # Step 4: Order the results by reward
    result_list.sort(key=lambda i: i["reward"])
    # Step 5: Add the n-last to the list of combinations to remove
    hit_list = result_list[:n_to_remove]
    for i in hit_list:
        print("Removing {0} with reward {1}".format(i["combination"], i["reward"]))
    # Step 6: remove the bad guys from the dataframe
    new_dataframe = dataframe_removal(dataframe, hit_list)
    return new_dataframe
