"""
    Fktable Loader using pineappl
"""
import numpy as np
import pandas as pd
from reportengine.compat import yaml
from validphys.coredata import FKTableData

EXT = "pineappl.lz4"


class PineAPPLEquivalentNotKnown(Exception):
    pass


class YamlFileNotFound(FileNotFoundError):
    pass


def _load_yaml(yaml_file):
    """Load a dataset.yaml file"""
    if not yaml_file.exists():
        raise YamlFileNotFound(yaml_file)
    ret = yaml.safe_load(yaml_file.read_text())
    # Make sure the operations are upper-cased for compound-compatibility
    ret["operation"] = "NULL" if ret["operation"] is None else ret["operation"].upper()
    return ret


def get_yaml_information(yaml_file, theorypath, check_pineappl=False):
    """
    Given a yaml_file, return the corresponding dictionary
    with all information and an extra field "paths"
    with all the grids to be loaded for the given dataset.
    Checks whether the grid is apfelcomb or pineappl
    if check_pineappl is True this function will raise PineAPPLEquivalentNotKnown
    if a pineappl grid is not found.
    Parameters
    ----------
        yaml_file: Path
            path of the yaml file for the given dataset
        theorypath: Path
            path of the theory folder where to find the grids
    Returns
    -------
        yaml_content: dict
            Metadata prepared for the FKTables
        paths: list(list(path))
            List (of lists) with all the grids that will need to be loaded
    """
    yaml_content = _load_yaml(yaml_file)
    grids_folder = theorypath / "pineappls"

    if yaml_content.get("appl") and check_pineappl:
        # This might be useful to use the "legacy loader" when there is no actual pineappl available
        raise PineAPPLEquivalentNotKnown(yaml_content["target_dataset"])

    # Turn the operans and the members into paths (and check all of them exist)
    ret = []
    for operand in yaml_content["operands"]:
        tmp = []
        for member in operand:
            p = grids_folder / f"{member}.{EXT}"
            if not p.exists():
                raise FileNotFoundError(f"Failed to find {p}")
            tmp.append(p)
        ret.append(tmp)

    return yaml_content, ret


def _pinelumi_to_vplumi(pine_luminosity, hadronic):
    """Convert the pineappl luminosity into vp-luminosity
    i.e., the non-zero indices of the 14x14 luminosity matrix for hadronic
    and the non-zero indices of the 14 flavours for DIS
    """
    eko_numbering_scheme = (22, 100, 21, 200, 203, 208, 215, 224, 235, 103, 108, 115, 124, 135)
    if hadronic:
        flavour_map = np.zeros((14, 14), dtype=bool)
        for i, j in pine_luminosity:
            idx = eko_numbering_scheme.index(i)
            jdx = eko_numbering_scheme.index(j)
            flavour_map[idx, jdx] = True
            co = np.where(flavour_map.ravel())[0]
    else:
        # The proton might come from both sides
        try:
            co = [eko_numbering_scheme.index(i) for _, i in pine_luminosity]
        except ValueError:
            co = [eko_numbering_scheme.index(i) for i, _ in pine_luminosity]
    return co


def pineappl_reader(fkspec):
    """
    Receives a fkspec, which contains the appropiate references to
    the pineappl grid to load and concatenate

    For more information: https://pineappl.readthedocs.io/en/latest/modules/pineappl/pineappl.html#pineappl.pineappl.PyFkTable
    """
    try:
        from pineappl.fk_table import FkTable
    except ImportError as e:
        raise ImportError("Please install pineappl: ~$ pip install pineappl") from e

    # Read each of the pineappl fktables
    pines = [FkTable.read(i) for i in fkspec.fkpath]

    # Extract some theory metadata from the first grid
    pine_representative = pines[0]
    hadronic = (
        pine_representative.key_values()["initial_state_1"]
        == pine_representative.key_values()["initial_state_2"]
    )
    Q0 = np.sqrt(pine_representative.muf2())
    xgrid = pine_representative.x_grid()
    # fktables in pineapplgrid are for o = fk * f while previous fktables were o = fk * xf
    # prepare the grid all tables will be divided by
    if hadronic:
        xdivision = (xgrid[:, None] * xgrid[None, :]).flatten()
    else:
        xdivision = xgrid

    ####
    # Prepare the apfelcomb-pineappl compatibility fixes
    # Note that these are applied at the grid level
    # (so before the fktable -concatenation of all grids- is constructed)
    operands = fkspec.metadata["operands"]
    apfelcomb_norm = None
    if fkspec.metadata.get("apfelcomb_norm"):
        norms = fkspec.metadata.get("apfelcomb_norm")  # would like to use the Walrus
        # There's no easy way for an fktable to know its role in given operation:
        for factors, grids in zip(norms, operands):
            # Now check, in case of an operation, what is our index in this operation
            if len(grids) == len(fkspec.fkpath) and all(
                f.name == f"{o}.{EXT}" for f, o in zip(fkspec.fkpath, grids)
            ):
                apfelcomb_norm = np.array(factors)

    # Check for the repetition flag, meaning we only want the first datapoint for this fktable
    apfelcomb_repetition_flag = False
    if fkspec.metadata.get("repetition_flag"):
        valid_targets = []
        for operand, flagged in zip(operands, fkspec.metadata["repetition_flag"]):
            if flagged:
                valid_targets.append(f"{operand[0]}.{EXT}")
        # Now check whether the current fktable is part of the valid targets
        apfelcomb_repetition_flag = fkspec.fkpath[0].name in valid_targets
        if apfelcomb_repetition_flag and len(fkspec.fkpath) > 1:
            raise ValueError(f"Repetition set for a group of fktables at once: {fkspec.fkpath}")

    # afaik there's only dataset with shifts, but it needs to be considered
    shifts = None
    if fkspec.metadata.get("shifts"):
        if len(operands) > 1:
            raise ValueError("Wrong shifts for {fkspec.metadata['target_dataset']}")
        shifts = [0 if shift is None else shift for shift in fkspec.metadata["shifts"][0]]
    ###

    # Read each separated grid and luminosity
    fktables = []
    ndata = 0
    for i, p in enumerate(pines):
        luminosity_columns = _pinelumi_to_vplumi(p.lumi(), hadronic)

        # Check whether this is a normalization table
        if fkspec.norm:
            raw_fktable = np.sum(p.table(), axis=0, keepdims=True)
        else:  # remove the bin normalization
            raw_fktable = (p.table().T / p.bin_normalizations()).T
        n = raw_fktable.shape[0]
        lf = len(luminosity_columns)

        # Apply the apfelcomb fixex
        if apfelcomb_norm is not None:
            raw_fktable = (raw_fktable.T * apfelcomb_norm[ndata : ndata + n]).T
        if apfelcomb_repetition_flag:
            raw_fktable = raw_fktable[0:1]
            n = 1
        if shifts is not None:
            ndata += shifts[i]
        ###

        partial_fktable = raw_fktable.reshape(n, lf, -1) / xdivision

        # Now concatenate (data, x1, x2) and move the flavours to the columns
        df_fktable = partial_fktable.swapaxes(0, 1).reshape(lf, -1).T

        # Create the multi-index for the dataframe
        ni = np.arange(ndata, n + ndata)
        xi = np.arange(len(xgrid))
        if hadronic:
            idx = pd.MultiIndex.from_product([ni, xi, xi], names=["data", "x1", "x2"])
        else:
            idx = pd.MultiIndex.from_product([ni, xi], names=["data", "x"])

        df_fktable *= fkspec.metadata.get("conversion_factor", 1.0)

        fktables.append(pd.DataFrame(df_fktable, columns=luminosity_columns, index=idx))
        ndata += n

    # Finallly concatenate all fktables, sort by flavours and fill any holes
    df = pd.concat(fktables, sort=True, copy=False).fillna(0.0)

    return FKTableData(
        sigma=df, ndata=ndata, Q0=Q0, metadata=fkspec.metadata, hadronic=hadronic, xgrid=xgrid
    )
