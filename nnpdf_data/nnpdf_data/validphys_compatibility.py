"""
    This file exists solely for me to be able to upload a package to PyPI for the nnpdf data which does not depend
    on the rest of the NNPDF code.
    This file should not be modified. Everything in this file is deprecated and should be removed, and the only reason
    this is needed is because we are still mixing new and old data.

    This also means that _things_ that would've been loaded as a kinematic transformation or result transformation
    are loaded as boring strings that cannot do anything, as it should be.
"""

import dataclasses
from functools import cache
import pathlib
import typing

import ruamel.yaml as yaml

labeler_functions = []


@dataclasses.dataclass
class PlottingOptions:
    func_labels: dict = dataclasses.field(default_factory=dict)
    dataset_label: typing.Optional[str] = None
    experiment: typing.Optional[str] = None
    nnpdf31_process: typing.Optional[str] = None
    data_reference: typing.Optional[str] = None
    theory_reference: typing.Optional[str] = None
    process_description: typing.Optional[str] = None
    y_label: typing.Optional[str] = None
    x_label: typing.Optional[str] = None
    kinematics_override: typing.Optional[str] = None
    result_transform: typing.Optional[str] = None
    x: typing.Optional[str] = None
    plot_x: typing.Optional[str] = None
    x_scale: typing.Optional[str] = None
    y_scale: typing.Optional[str] = None
    line_by: typing.Optional[list] = None
    figure_by: typing.Optional[list] = None
    extra_labels: typing.Optional[typing.Mapping[str, typing.List]] = None
    normalize: typing.Optional[dict] = None
    # Note that this "PlottingOptions" start already digested, because it actually does nothing!
    already_digested: typing.Optional[bool] = True


########## Legacy names compatibility
# The functions and variables below exist solely for reproducibility of old NNPDF results
# and no external code or new feature should depend on them as they might be removed at any point
# with no previous warning

path_vpdata = pathlib.Path(__file__).parent
path_commondata = path_vpdata / "commondata"
_path_legacy_mapping = path_commondata / "dataset_names.yml"
_legacy_to_new_mapping_raw = yaml.YAML().load(_path_legacy_mapping)
# Convert strings into a dictionary
legacy_to_new_mapping = {
    k: ({"dataset": v} if isinstance(v, str) else v) for k, v in _legacy_to_new_mapping_raw.items()
}


@cache
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


@cache
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
