#!/usr/bin/env python
"""
Script to convert old commondata into the new format

Note that the uncertainties for autoconverted sets are always set as legacy
"""
from argparse import ArgumentParser
import functools
from pathlib import Path
import traceback

import pandas as pd
from yaml import safe_dump, safe_load

# Prepare the folders by automatically finding them from validphys
from validphys.datafiles import path_commondata as old_cd_root
from validphys.loader import Loader

# old_cd_root = validphys.datafiles / "commondata"
new_cd_root = old_cd_root.with_name("new_commondata")
dataset_names_path = new_cd_root / "dataset_names.yml"
theory_files = Loader()._theories_path


def read_commondata_csv(commondatafile):
    """Read the old format commondata which were csv files we sparkling formatting
    This is directly taken from validphys
    """
    commondatatable = pd.read_csv(commondatafile, sep=r"\s+", skiprows=1, header=None)
    # Do we have NaNs? files with wrong formatting?
    commondataheader = ["entry", "process", "kin1", "kin2", "kin3", "data", "stat"]
    nsys = (commondatatable.shape[1] - len(commondataheader)) // 2

    commondataheader += ["ADD", "MULT"] * nsys
    commondatatable.columns = commondataheader
    commondatatable.set_index("entry", inplace=True)
    return commondatatable


def create_data(commondata_df):
    """Given a commondata dataframe, extract the central data and create a dictionary
    ready to be yamld'd"""
    data = commondata_df["data"].values
    return {"data_central": data.tolist()}


def create_kinematics(df):
    """Create kinematics dictionary with kin1, kin2, kin3 . . ."""
    kin_df = df[["kin1", "kin2", "kin3"]]
    bins = []
    for _, b in kin_df.T.items():
        tmp = {}
        for k, val in b.items():
            tmp[f"k{k[-1]}"] = {"min": None, "mid": val, "max": None}
        bins.append(tmp)
    return {"bins": bins}


def create_uncertainties(df, systype_file, is_default=False, use_multiplicative=False):
    """Create the uncertainties dictionary from the old cd information

    We first clean the dataframe to have only the systematic uncertainties
    and then even-indexes will be ADD uncertainties and odd-indexes MULT uncertainties

    Such that e.g., unc 4 in the systype file, if it is MULT, will correspond to index 7
    """
    stat = df["stat"].values.tolist()
    data = df["data"].values

    to_drop = ["process", "kin1", "kin2", "kin3", "data", "stat"]
    unc_df = df.drop(to_drop, axis=1)

    sys_df = pd.read_csv(systype_file, sep=r"\s+", skiprows=1, header=None, index_col=0)
    definitions = {}

    for i, unc_type in sys_df.T.items():
        definitions[f"sys_corr_{i}"] = {
            "description": f"Sys uncertainty idx: {i}",
            "treatment": unc_type[1],
            "type": unc_type[2],
        }

    # Check whether the number of uncertainties in the systype file is consistent with the df
    if len(unc_df.columns) != len(sys_df) * 2:
        # If this happened and this is the (true) DEFAULT, crash
        if is_default:
            raise ValueError("Different number of systematics in systype file and commondata")
        return {}

    bins = []
    for idx, bin_data in unc_df.T.items():
        tmp = {"stat": stat[idx - 1]}
        for n, (key, info) in enumerate(definitions.items()):
            if info["treatment"] not in ["ADD", "MULT"]:
                raise ValueError(f"Treatment type: {info['treatment']} not recognized")
            if use_multiplicative:
                tmp[key] = float(bin_data[2 * n + 1] * data[idx - 1] / 100.0)
            else:
                tmp[key] = float(bin_data[2 * n])
        bins.append(tmp)

    # Now add stat to the definitions
    definitions_out = {
        "stat": {
            "description": "Uncorrelated statistical uncertainties",
            "treatment": "ADD",
            "type": "UNCORR",
        },
        **definitions,
    }

    return {"definitions": definitions_out, "bins": bins}


def create_plotting(plotting_file, plotting_type=None):
    """Create the plotting dictionary by merging the plotting file of the commondata
    and the plotting type associated to it
    """
    type_info = {}
    if plotting_type is not None and plotting_type.exists():
        type_info = safe_load(plotting_type.read_text())

    plotting_data = safe_load(plotting_file.read_text())

    plotting_dict = {**type_info, **plotting_data}

    plotting_dict["plot_x"] = plotting_dict.pop("x", 'idat')
    return plotting_dict


def create_theory(yamldb_files):
    """Create the theory dictionary by reading the yamldb file (if available)"""
    if not yamldb_files:
        return {}

    # Assume the last one is the best most updated one:
    yaml_info = safe_load(yamldb_files[-1].read_text())
    # Change "operands" into "fktables" and drop target
    yaml_info.pop("target_dataset", None)
    yaml_info.pop("comment", None)
    yaml_info["FK_tables"] = yaml_info.pop("operands")

    operation = yaml_info.get("operation")
    if operation is None:
        operation = 'NULL'
    yaml_info["operation"] = operation

    return yaml_info


def _gen_k(kx):
    """Generate k{x} variable"""
    return {"description": f"Variable {kx}", "label": f"{kx}", "units": ""}


def create_obs_dict(commondata_df, plotting_dict, theory_dict, obs_name="PLACEHOLDER"):
    """Create the observable dictionary by combining available information
    in the commondata dataframe and the plotting_dict

    It doesn't fill any data files (i.e., uncertainties, data or kinematics)
    """
    final_plotting_dict = dict(plotting_dict)

    # Extract necessary information
    ndata = len(commondata_df)
    process_type = commondata_df["process"][1]

    description = final_plotting_dict.pop("process_description", "DESCRIPTION_PLACEHOLDER")
    label = final_plotting_dict["dataset_label"]
    units = ""

    final_plotting_dict.pop("nnpdf31_process")
    final_plotting_dict.pop("experiment")

    # Sub-dicts
    observable = {"description": description, "label": label, "units": units}

    coverage = ["k1", "k2", "k3"]
    kinematics = {"variables": {i: _gen_k(i) for i in coverage}}

    return {
        "observable_name": obs_name,
        "observable": observable,
        "process_type": process_type,
        "tables": [],
        "npoints": [],
        "ndata": ndata,
        "plotting": final_plotting_dict,
        "kinematic_coverage": coverage,
        "kinematics": kinematics,
        "theory": theory_dict,
        "data_uncertainties": [],
    }


def yaml_dump_wrapper(data, target_file, dry=False, **kwargs):
    """Wrapper around safe_sump in order to use the dry flag"""
    if dry:
        return None
    safe_dump(data, target_file.open("w", encoding="utf-8"), **kwargs)


def get_all_data_files(dsname):
    """Given the (old) name of a dataset, list the data file,
    the systematics default file and possible other systematics"""
    data_file = old_cd_root / f"DATA_{dsname}.dat"
    if not data_file.exists():
        raise FileNotFoundError(f"commondata file {data_file} not found")

    # systypes
    systypes_default = old_cd_root / "systypes" / f"SYSTYPE_{dsname}_DEFAULT.dat"
    if not systypes_default.exists():
        raise FileNotFoundError(f"Not default systypes found")

    # Other systypes?
    other_systypes = list((old_cd_root / "systypes").glob(f"SYSTYPE_{dsname}_*"))
    # Remove the default
    other_systypes.remove(systypes_default)
    return data_file, systypes_default, other_systypes


def convert_from_old_to_new(dsname, new_info, overwrite=False, dry=False, keep_existing=False):
    """
    Convert the old commondata ``dsname``
    into the new commondata format defined by ``new_info`` which includes:
        ``{new_name, reference}``

    If an entry already exists for the old commondata in the ``dataset_names.yaml`` file
    then exit without doing anything.
    For now we assume a manually written dataset trumps any automatic implementation.
    If the overwrite flag is set to True, overwrite only if the correspondance between
    old and new is through a ``legacy`` variant.

    If an entry for the new commondata already exist in the target metadata file
    then raise an Exception unless the overwrite flag is explicitly set to True.

    Positivity datasets (defined because their experiment name is NNPDF and their process POS)
    are treated slightly different

    If overwrite is False, the script will raise an Exception if a metadata
    already exists and any inconsistency is found between it and the new data

    If overwrite is true and keep_existing is also true, when an observable already exit,
    extra information will be added to the legacy dictionary, but nothing will be overwritten
    """
    yaml_safe_dump = functools.partial(yaml_dump_wrapper, dry=dry)

    new_name = new_info["new_name"]
    reference_arxiv = new_info.get("reference_arxiv", "")
    reference_journal = new_info.get("journal", "")
    reference_hepdata = new_info.get("reference_hepdata", "")

    variant_name = new_info.get("variant_name", "legacy")

    # Is this an special positivity dataset?
    is_positivity_ds = new_name.startswith(("NNPDF_POS", "NNPDF_INTEG"))

    dataset_info = safe_load(dataset_names_path.read_text())
    if dsname in dataset_info:
        if not (dataset_info[dsname].get("variant") == variant_name and overwrite):
            print(f"An entry for {dsname} already exist, skipping")
            return 0

    data_file, systypes_default, other_systypes = get_all_data_files(dsname)
    # If systypes exist, go through them to create extra (legacy, or another name,) variants
    extra_variants = []

    # Any yamldb files?
    yamldb_files = list(theory_files.glob(f"*/yamldb/{dsname}.yaml"))

    # plotting files
    plotting_file = old_cd_root / f"PLOTTING_{dsname}.yaml"
    if not plotting_file.exists():
        # Try the yml version
        if not plotting_file.with_suffix(".yml").exists():
            raise FileNotFoundError(f"No plotting file found: {plotting_file}")
        plotting_file = plotting_file.with_suffix(".yml")

    commondata_df = read_commondata_csv(data_file)
    proc = commondata_df["process"][1]

    plotting_type = old_cd_root / f"PLOTTINGTYPE_{proc}.yaml"
    # if the plotting type is DIS_CC or DIS_NC, the plotting_type is DIS
    if proc[:3] == "DIS":
        # Same check that it is done in validphys
        plotting_type = old_cd_root / f"PLOTTINGTYPE_DIS.yaml"
    if proc[:3] == "DYP":
        # Same check that it is done in validphys
        plotting_type = old_cd_root / f"PLOTTINGTYPE_DYP.yaml"
    if proc == "EWK_RAP_ASY":
        # oh my
        plotting_type = old_cd_root / f"PLOTTINGTYPE_EWK_RAP.yaml"

    # Now create the information that will be saved into the new commondata yaml files
    kinematics_dict = create_kinematics(commondata_df)
    plotting_dict = create_plotting(plotting_file, plotting_type)
    theory_dict = create_theory(yamldb_files)

    # Here the logic changes for positivity datasets
    if is_positivity_ds:
        uncertainties_dict = {}
        plotting_dict["experiment"] = "NNPDF"
        if new_name.startswith("NNPDF_INTEG"):
            plotting_dict["nnpdf31_process"] = "POSITIVITY"
        else:
            plotting_dict["nnpdf31_process"] = "INTEGRABILITY"
    else:
        data_dict = create_data(commondata_df)
        uncertainties_dict = create_uncertainties(commondata_df, systypes_default, is_default=True)

        # Note that some systypes might be "fake" variants due to the fact that some are also data variants
        # that granularity will be solved by hand unless the number of uncertainties is different
        for systype_file in other_systypes:
            tmp = create_uncertainties(commondata_df, systype_file)
            if tmp:
                extra_variant_name = systype_file.name.replace(f"SYSTYPE_{dsname}_", "").replace(
                    ".dat", ""
                )
                extra_variants.append((extra_variant_name, tmp))

    # Separate set name and observable
    if (set_name := new_info.get("set_name")) is None:
        obs_name = new_name.rsplit("_", 1)[-1]
        set_name = new_name.replace(f"_{obs_name}", "")
    else:
        obs_name = new_name.replace(f"{set_name}_", "")

    # Create the observable dictionary for this dataset
    obs_dict = create_obs_dict(commondata_df, plotting_dict, theory_dict, obs_name=obs_name)

    # Create the folder where to put the files
    output_folder = new_cd_root / set_name
    output_folder.mkdir(exist_ok=True, parents=True)

    # Read metadata (if it exist!) and perform necessary sanity checks before writing anything down
    metadata_path = output_folder / "metadata.yaml"
    existing_dataset = None

    if metadata_path.exists():
        metadata = safe_load(metadata_path.read_text())
        # Perform sanity checks
        nnpdf_md = metadata["nnpdf_metadata"]
        try:
            assert nnpdf_md["experiment"] == plotting_dict["experiment"]
            assert nnpdf_md["nnpdf31_process"] == plotting_dict["nnpdf31_process"]
            assert metadata.get("setname") == set_name
        except AssertionError:
            print(traceback.format_exc())
            # If this fails, inspect
            import ipdb

            ipdb.set_trace()
        # Check whether the observable already exists
        already_implemented = [i["observable_name"] for i in metadata["implemented_observables"]]
        if obs_name in already_implemented:
            if not overwrite:
                raise ValueError(f"{obs_name} already implemented for {set_name}")
            idx = already_implemented.index(obs_name)

            existing_dataset = metadata["implemented_observables"].pop(idx)
            if not keep_existing:
                existing_dataset = None
    else:
        # Create it anew!
        nnpdf_md = {
            "nnpdf31_process": plotting_dict["nnpdf31_process"],
            "experiment": plotting_dict["experiment"],
        }
        metadata = {
            "setname": set_name,
            "version": 1,
            "version_comment": "Port of old commondata",
            "nnpdf_metadata": nnpdf_md,
            "arXiv": {"url": f"https://arxiv.org/abs/{reference_arxiv}"},
            "iNSPIRE": {"url": ""},
            "hepdata": {"url": reference_hepdata, "version": -1},
            "implemented_observables": [],
        }
        if reference_journal is not None and reference_journal.strip():
            metadata["arXiv"]["journal"] = reference_journal

    # Are we keeping information, if so, we need a sanity check...
    if existing_dataset is not None:
        # It is truly the same data?
        # Needs to have the same number of datapoints!
        if existing_dataset["ndata"] != obs_dict["ndata"]:
            raise ValueError(
                "The dataset {old_name} has a different number of datapoints that {new_name}"
            )

        # Note that by doing this we are assuming that the theory already set is the best one
        # most of the time it will be correct
        obs_dict = existing_dataset
    else:
        # Put the files in the folder and update the observable dictionary
        kin_path = output_folder / f"kinematics_{obs_name}.yaml"
        yaml_safe_dump(kinematics_dict, kin_path, sort_keys=False)
        obs_dict["kinematics"]["file"] = kin_path.name

    dataset_info[dsname] = {"dataset": new_name}

    if not is_positivity_ds:
        data_path = output_folder / f"data_{variant_name}_{obs_name}.yaml"
        unc_path = output_folder / f"uncertainties_{variant_name}_{obs_name}.yaml"

        dataset_info[dsname]["variant"] = variant_name

        yaml_safe_dump(data_dict, data_path)
        yaml_safe_dump(uncertainties_dict, unc_path, sort_keys=False)

        # Put the uncertainties under "variants"
        if "variants" in obs_dict:
            obs_dict["variants"][variant_name] = {"data_uncertainties": [unc_path.name]}
        else:
            obs_dict["variants"] = {variant_name: {"data_uncertainties": [unc_path.name]}}

        # If this is a new dataset, use `data_central`, otherwise put it under variants
        if existing_dataset is None:
            obs_dict["data_central"] = data_path.name
        else:
            obs_dict["variants"][variant_name]["data_central"] = data_path.name

    # Add an extra key
    obs_dict["ported_from"] = dsname

    # Now loop over possible extra variants for different systypes
    for extra_variant_name, variant_dict in extra_variants:
        var_path = output_folder / f"uncertainties_{obs_name}_sys_{extra_variant_name}.yaml"
        yaml_safe_dump(variant_dict, var_path, sort_keys=False)
        obs_dict["variants"][f"legacy_{extra_variant_name}"] = {
            "data_uncertainties": [var_path.name]
        }

    metadata["implemented_observables"].append(obs_dict)
    yaml_safe_dump(metadata, metadata_path, sort_keys=False)

    # Update dataset_names
    yaml_safe_dump(dataset_info, dataset_names_path, sort_keys=False)

    print(f"Written new cd for {set_name}_{obs_name} to {output_folder}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("mapping_file", help="Mapping file with old-new information", type=Path)
    parser.add_argument(
        "--overwrite", help="Overwrite previous existing observable", action="store_true"
    )
    parser.add_argument(
        "--keep",
        help="To be used together with overwrite, whether to keep existing information",
        action="store_true",
    )
    parser.add_argument("--dry", help="Don't change any files", action="store_true")
    parser.add_argument("--filter", help="Apply a filter on the old name", type=str)
    args = parser.parse_args()

    mapping_info = safe_load(args.mapping_file.read_text())
    for old_name, new_info in mapping_info.items():
        if args.filter is not None and args.filter not in old_name:
            continue
        convert_from_old_to_new(
            old_name, new_info, overwrite=args.overwrite, dry=args.dry, keep_existing=args.keep
        )
