"""
    Loader for the pineappl-based FKTables

    The FKTables for pineappl have ``pineappl.lz4`` and can be utilized
    directly with the ``pineappl`` cli as well as read with ``pineappl.fk_table``
"""

import logging

import numpy as np
import pandas as pd

from validphys.commondataparser import EXT, TheoryMeta
from validphys.coredata import FKTableData

log = logging.getLogger(__name__)


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
    theory_meta = TheoryMeta.parser(yaml_file)
    member_paths = theory_meta.fktables_to_paths(grids_folder)
    return theory_meta, member_paths


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

    pines = []
    for fk_path in fkspec.fkpath:
        try:
            pines.append(FkTable.read(fk_path))
        except BaseException as e:
            # Catch absolutely any error coming from pineappl, give some info and immediately raise
            log.error(f"Fatal error reading {fk_path}")
            raise e

    cfactors = fkspec.load_cfactors()

    # Extract metadata from the first grid
    pine_rep = pines[0]

    # Check if it is Polarized FK table
    is_polarized = pine_rep.key_values().get("polarized") == "True"

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

    # Process the shifts and normalizations (if any),
    # shifts is a dictionary with {fktable_name: shift_value}
    # normalization instead {fktable_name: normalization to apply}
    # since this parser doesn't know about operations, we need to convert it to a list
    # then we just iterate over the fktables and apply the shift in the right order
    shifts = fkspec.metadata.shifts
    normalization_per_fktable = fkspec.metadata.normalization
    fknames = [i.name.replace(f".{EXT}", "") for i in fkspec.fkpath]
    if cfactors is not None:
        cfactors = dict(zip(fknames, cfactors))

    # fktables in pineapplgrid are for obs = fk * f while previous fktables were obs = fk * xf
    # prepare the grid all tables will be divided by
    if hadronic:
        xdivision = np.prod(np.meshgrid(xgrid, xgrid), axis=0)
    else:
        xdivision = xgrid[:, np.newaxis]

    partial_fktables = []
    ndata = 0
    for fkname, p in zip(fknames, pines):
        # Start by reading possible cfactors if cfactor is not empty
        cfprod = 1.0
        if cfactors is not None:
            for cfac in cfactors.get(fkname, []):
                cfprod *= cfac.central_value

        # Read the table, remove bin normalization and apply cfactors
        raw_fktable = (cfprod * p.table().T / p.bin_normalizations()).T
        n = raw_fktable.shape[0]

        # Apply possible per-fktable fixes
        if shifts is not None:
            ndata += shifts.get(fkname, 0)

        if normalization_per_fktable is not None:
            raw_fktable = raw_fktable * normalization_per_fktable.get(fkname, 1.0)

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

    # Check whether this is a 1-point normalization fktable and, if that's the case, protect!
    if fkspec.metadata.operation == "RATIO" and len(pines) == 1:
        # it _might_ be, check whether it is the divisor fktable
        divisor = fkspec.metadata.FK_tables[-1][0]
        name = fkspec.fkpath[0].name.replace(f".{EXT}", "")

        if np.allclose(sigma.loc[1:], 0.0):
            # Old denominator fktables were filled with 0s beyond the first point
            # and they needed to be post-processed to repeat the same point many time
            # Instead, drop everything beyond the 1st point (used 0:0 to keep the same kind of df)
            sigma = sigma.loc[0:0]
            ndata = 1

        if ndata == 1:
            # There's no doubt
            protected = divisor == name

    return FKTableData(
        sigma=sigma,
        ndata=ndata,
        Q0=Q0,
        is_polarized=is_polarized,
        metadata=fkspec.metadata,
        hadronic=hadronic,
        xgrid=xgrid,
        protected=protected,
    )
