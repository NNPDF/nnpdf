from functools import lru_cache
import pathlib

import ruamel.yaml as yaml

from ._version import __version__

path_vpdata = pathlib.Path(__file__).parent
path_commondata = path_vpdata / "commondata"

# VP should not have access to this file, only to the products
_path_legacy_mapping = path_commondata / "dataset_names.yml"
theory_cards = path_vpdata / "theory_cards"

_legacy_to_new_mapping_raw = yaml.YAML().load(_path_legacy_mapping)
# Convert strings into a dictionary
legacy_to_new_mapping = {
    k: ({"dataset": v} if isinstance(v, str) else v) for k, v in _legacy_to_new_mapping_raw.items()
}


@lru_cache
def legacy_to_new_map(dataset_name, sys=None):
    """Find the new dataset name and variant corresponding to an old dataset
    and systematics choice"""
    if dataset_name not in legacy_to_new_mapping:
        return dataset_name, None

    new_name = legacy_to_new_mapping[dataset_name]
    variant = new_name.get("variant")
    new_name = new_name["dataset"]
    if sys is not None:
        if variant is None:
            raise KeyError(
                f"I cannot translate the combination of {dataset_name} and sys: {sys}. Please report this."
            )
        variant += f"_{sys}"

    return new_name, variant


@lru_cache
def new_to_legacy_map(dataset_name, variant_used):
    """Loop over the dictionary and find the right dataset.

    Since it is posible to have more than 1 dataset mapped to the same new one,
    returns a list of everything that matches.

    This function will loop over the entire dictionary of mappings and selects
    1. All datasets that match exactly what's in the runcard (dataset & variant): exact_matches
    2. All datasets that match the dataset name: matches
    If there are any `exact_matches`, it will return only those; otherwise, return all `matches`
    if there are no `matches` at all, return None
    """

    matches = []
    exact_matches = []

    for old_name, new_info in legacy_to_new_mapping.items():
        new_name = new_info["dataset"]
        variant = new_info.get("variant")

        if new_name == dataset_name:
            matches.append(old_name)
            # if it's a nuclear DIS data promote legacy to be legacy_dw
            if "_DW_" in old_name and variant_used == "legacy":
                variant = "legacy_dw"

            if variant_used == variant:
                exact_matches.append(old_name)

    # If we found exact matches, return those and stop looking
    if exact_matches:
        return exact_matches
    elif matches:
        return matches
    return None
