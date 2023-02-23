"""
    Loader for the pineappl-based FKTables

    The FKTables for pineappl have ``pineappl.lz4`` and can be utilized
    directly with the ``pineappl`` cli as well as read with ``pineappl.fk_table``
"""
import numpy as np
import pandas as pd

from reportengine.compat import yaml
from validphys.coredata import FKTableData
from validphys.commondataparser import TheoryMeta, EXT
from validphys.utils import parse_yaml_inp


class GridFileNotFound(FileNotFoundError):
    """PineAPPL file for FK table not found."""


def pineko_yaml(yaml_file, grids_folder):
    """Given a yaml_file, returns the corresponding dictionary and grids.

    The dictionary contains all information and we return an extra field
    with all the grids to be loaded for the given dataset.

    Parameters
    ----------
    yaml_file : pathlib.Path
        path of the yaml file for the given dataset
    grids_folder : pathlib.Path
        path of the grids folder
    check_grid_existence: bool
        if True (default) checks whether the grid exists

    Returns
    -------
    yaml_content: dict
        Metadata prepared for the FKTables
    paths: list(list(path))
        List (of lists) with all the grids that will need to be loaded
    """
    # TODO: the theory metadata can be found inside the commondata metadata
    # however, for the time being, pineappl tables contain this information in the `yamldb` database
    # they should be 100% compatible (and if they are not there is something wrong somewhere)
    # so already at this stage, use TheoryMeta parser to get the metadata for pineappl theories
    # Note also that we need to use this "parser" due to the usage of the name "operands" in the yamldb
#     theory_meta = parse_yaml_inp(yaml_file, TheoryMeta)
    theory_meta = TheoryMeta.parser(yaml_file)
    member_paths = theory_meta.fktables_to_paths(grids_folder)
    return theory_meta, member_paths


def pineko_apfelcomb_compatibility_flags(gridpaths, metadata):
    """
    Prepare the apfelcomb-pineappl compatibility fixes by matching the apfelcomb content
    of the metadata to the grids that are being loaded.

    These fixes can be of only three types:

    - normalization:
        normalization per subgrid

        normalization:
            grid_name: factor

    - repetition_flag:
        when a grid was actually the same point repeated X times
        NNPDF cfactors and cuts are waiting for this repetition and so we need to keep track of it

        repetition_flag:
            grid_name

    - shifts:
        only for ATLASZPT8TEVMDIST
        the points in this dataset are not contiguous so the index is shifted

        shifts:
            grid_name: shift_int

    Returns
    -------
        apfelcomb_norm: np.array
            Per-point normalization factor to be applied to the grid
            to be compatible with the data
        apfelcomb_repetition_flag: bool
            Whether the fktable is a single point which gets repeated up to a certain size
            (for instance to normalize a distribution)
        shift: list(int)
            Shift in the data index for each grid that forms the fktable
    """
    apfelcomb = metadata.apfelcomb
    if apfelcomb is None:
        return None

    # Can't pathlib understand double suffixes?
    operands = [i.name.replace(f".{EXT}", "") for i in gridpaths]
    ret = {}

    # Check whether we have a normalization active and whether it affects any of the grids
    if apfelcomb.normalization is not None:
        norm_info = apfelcomb.normalization
        # Now fill the operands that need normalization
        ret["normalization"] = [norm_info.get(op, 1.0) for op in operands]

    # Check whether the repetition flag is active
    if apfelcomb.repetition_flag is not None:
        if len(operands) == 1:
            ret["repetition_flag"] = operands[0] in apfelcomb.repetition_flag
        else:
            # Just for the sake of it, let's check whether we did something stupid
            if any(op in apfelcomb.repetition_flag for op in operands):
                raise ValueError(f"The yaml info for {metadata['target_dataset']} is broken")

    # Check whether the dataset has shifts
    # NOTE: this only happens for ATLASZPT8TEVMDIST, if that gets fixed we might as well remove it
    if apfelcomb.shifts is not None:
        shift_info = apfelcomb.shifts
        ret["shifts"] = [shift_info.get(op, 0) for op in operands]

    return ret


def _pinelumi_to_columns(pine_luminosity, hadronic):
    """Makes the pineappl luminosity into the column indices of a dataframe
    These corresponds to the indices of a flattened (14x14) matrix for hadronic observables
    and the non-zero indices of the 14-flavours for DIS

    Parameters
    ----------
        pine_luminosity: list(tuple(int))
            list with a pair of flavours per channel
        hadronic: bool
            flag for hadronic / DIS observables

    Returns
    -------
        list(int): list of labels for the columns
    """
    evol_basis_pids = tuple(
        [22, 100, 21, 200]
        + [200 + n**2 - 1 for n in range(2, 6 + 1)]
        + [100 + n**2 - 1 for n in range(2, 6 + 1)]
    )
    flav_size = len(evol_basis_pids)
    columns = []
    if hadronic:
        for i, j in pine_luminosity:
            idx = evol_basis_pids.index(i)
            jdx = evol_basis_pids.index(j)
            columns.append(flav_size * idx + jdx)
    else:
        # The proton might come from both sides
        try:
            columns = [evol_basis_pids.index(i) for _, i in pine_luminosity]
        except ValueError:
            columns = [evol_basis_pids.index(i) for i, _ in pine_luminosity]
    return columns


def get_yaml_information(yaml_file, theorypath):
    """Reads the yaml information from a yaml compound file

    Transitional function: the call to "pineko" might be to some other commondata reader
    that will know how to extract the information from the commondata
    """
    # The pineappl grids are just stored where the fktables would be, "fastkernel"
    grids_folder = theorypath / "fastkernel"
    return pineko_yaml(yaml_file, grids_folder)


def pineappl_reader(fkspec):
    """
    Receives a fkspec, which contains the path to the fktables that are to be read by pineappl
    as well as metadata that fixes things like conversion factors or apfelcomb flag.
    The fkspec contains also the cfactors which are applied _directly_ to each of the fktables.

    The output of this function is an instance of FKTableData which can be generated from reading
    several FKTable files which get concatenated on the ndata (bin) axis.

    For more information on the reading of pineappl tables:
        https://pineappl.readthedocs.io/en/latest/modules/pineappl/pineappl.html#pineappl.pineappl.PyFkTable

    About the reader:
        Each pineappl table is a 4-dimensional grid with:
            (ndata, active channels, x1, x2)
        for DIS grids x2 will contain one single number.
        The luminosity channels are given in a (flav1, flav2) format and thus need to be converted
        to the 1-D index of a (14x14) luminosity tensor in order to put in the form of a dataframe.

        All grids in pineappl are constructed with the exact same xgrid,
        the active channels can vary and so when grids are concatenated for an observable
        the gaps are filled with 0s.

        The pineappl grids are such that obs = sum_{bins} fk * f (*f) * bin_w
        so in order to use them together with old-style grids (obs = sum_{bins} fk * xf (*xf))
        it is necessary to remove the factor of x and the normalization of the bins.

    About apfelcomb flags in yamldb files:
        old commondata files and old grids have over time been through various iterations while remaining compatibility between each other,
        and fixes and hacks have been incorporated in one or another
        for the new theory to be compatible with old commpondata it is necessary
        to keep track of said hacks (and to apply conversion factors when required)
    NOTE: both conversion factors and apfelcomb flags will be eventually removed.

    Returns
    -------
        validphys.coredata.FKTableData
            an FKTableData object containing all necessary information to compute predictions
    """
    from pineappl.fk_table import FkTable

    pines = [FkTable.read(i) for i in fkspec.fkpath]
    cfactors = fkspec.load_cfactors()

    # Extract metadata from the first grid
    pine_rep = pines[0]

    # Is it hadronic? (at the moment only hadronic and DIS are considered)
    hadronic = pine_rep.key_values()["initial_state_1"] == pine_rep.key_values()["initial_state_2"]
    # Sanity check (in case at some point we start fitting things that are not protons)
    if hadronic and pine_rep.key_values()["initial_state_1"] != "2212":
        raise ValueError(
            "pineappl_reader is not prepared to read a hadronic fktable with no protons!"
        )
    Q0 = np.sqrt(pine_rep.muf2())
    xgrid = np.array([])
    for pine in pines:
        xgrid = np.union1d(xgrid, pine.x_grid())
    xi = np.arange(len(xgrid))
    protected = False

    apfelcomb = pineko_apfelcomb_compatibility_flags(fkspec.fkpath, fkspec.metadata)

    # fktables in pineapplgrid are for obs = fk * f while previous fktables were obs = fk * xf
    # prepare the grid all tables will be divided by
    if hadronic:
        xdivision = np.prod(np.meshgrid(xgrid, xgrid), axis=0)
    else:
        xdivision = xgrid[:, np.newaxis]

    partial_fktables = []
    ndata = 0
    for i, p in enumerate(pines):
        # Start by reading possible cfactors if cfactor is not empty
        cfprod = 1.0
        if cfactors:
            for cfac in cfactors[i]:
                cfprod *= cfac.central_value

        # Read the table, remove bin normalization and apply cfactors
        raw_fktable = (cfprod * p.table().T / p.bin_normalizations()).T
        n = raw_fktable.shape[0]

        # Apply the apfelcomb fixes _if_ they are needed
        if apfelcomb is not None:
            if apfelcomb.get("normalization") is not None:
                raw_fktable = raw_fktable * apfelcomb["normalization"][i]
            if apfelcomb.get("repetition_flag", False):
                raw_fktable = raw_fktable[0:1]
                n = 1
                protected = True
            if apfelcomb.get("shifts") is not None:
                ndata += apfelcomb["shifts"][i]

        # Add empty points to ensure that all fktables share the same x-grid upon convolution
        missing_x_points = np.setdiff1d(xgrid, p.x_grid(), assume_unique=True)
        for x_point in missing_x_points:
            miss_index = list(xgrid).index(x_point)
            raw_fktable = np.insert(raw_fktable, miss_index, 0.0, axis=2)
            if hadronic:
                raw_fktable = np.insert(raw_fktable, miss_index, 0.0, axis=3)
        # Check conversion factors and remove the x* from the fktable
        raw_fktable *= fkspec.metadata.conversion_factor / xdivision

        # Create the multi-index for the dataframe
        # for optimized pineappls different grids can potentially have different indices
        # so they need to be indexed separately and then concatenated only at the end
        lumi_columns = _pinelumi_to_columns(p.lumi(), hadronic)
        lf = len(lumi_columns)
        data_idx = np.arange(ndata, ndata + n)
        if hadronic:
            idx = pd.MultiIndex.from_product([data_idx, xi, xi], names=["data", "x1", "x2"])
        else:
            idx = pd.MultiIndex.from_product([data_idx, xi], names=["data", "x"])

        # Now concatenate (data, x1, x2) and move the flavours to the columns
        df_fktable = raw_fktable.swapaxes(0, 1).reshape(lf, -1).T
        partial_fktables.append(pd.DataFrame(df_fktable, columns=lumi_columns, index=idx))

        ndata += n

    # Finallly concatenate all fktables, sort by flavours and fill any holes
    sigma = pd.concat(partial_fktables, sort=True, copy=False).fillna(0.0)

    return FKTableData(
        sigma=sigma,
        ndata=ndata,
        Q0=Q0,
        metadata=fkspec.metadata,
        hadronic=hadronic,
        xgrid=xgrid,
        protected=protected,
    )
