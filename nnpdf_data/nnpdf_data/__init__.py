from functools import lru_cache
import pathlib

import ruamel.yaml as yaml

from ._version import __version__

path_vpdata = pathlib.Path(__file__).parent
path_commondata = path_vpdata / "new_commondata"

# VP should not have access to this file, only to the products
_path_legacy_mapping = path_commondata / "dataset_names.yml"
legacy_to_new_mapping = yaml.YAML().load(_path_legacy_mapping)
theory_cards = path_vpdata / "theory_cards"


@lru_cache
def legacy_to_new_map(dataset_name, sys=None):
    """Find the new dataset name and variant corresponding to an old dataset
    and systematics choice"""
    if dataset_name not in legacy_to_new_mapping:
        return dataset_name, None

    new_name = legacy_to_new_mapping[dataset_name]
    if isinstance(new_name, str):
        if sys is not None:
            raise KeyError(
                f"I cannot translate the combination of {dataset_name} and sys: {sys}. Please report this."
            )
        return new_name, None

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
    """Loop over the dictionary and find the right dataset"""
    # It is not possible to reverse the dictionary because
    # we can have 2 old dataset mapped to the same new one

    possible_match = None
    for old_name, new_name in legacy_to_new_mapping.items():
        variant = None
        if not isinstance(new_name, str):
            variant = new_name.get("variant")
            new_name = new_name["dataset"]

        if new_name == dataset_name:
            if variant_used == variant:
                return old_name
            # Now, for legacy variants we might want to match (sys,)
            # so accept anything that starts with `legacy_`
            # so variant `legacy_10` will match `legacy` in the dictionary
            # but if an exact match if found before, the search ends
            if variant_used is not None and variant_used.startswith("legacy_"):
                possible_match = old_name

    return possible_match
