"""
n3fit_data_utils.py

Library of function that read validphys object and make them into clean numpy objects
for its usage for the fitting framework
"""
# TODO:
# Wrap them into a single dataclass with all relevant information instead of using dictionaries
import numpy as np
import pandas as pd


def fk_parser(fk, is_hadronic=False):
    """
    # Arguments:
        - `fk`: fktable object

    # Return:
        - `dict_out`: dictionary with all information about the fktable
            - 'xgrid'
            - 'nx'
            - 'ndata'
            - 'basis'
            - 'fktable'
    """

    # Dimensional information
    nonzero = fk.GetNonZero()
    ndata = fk.GetNData()
    nx = fk.GetTx()

    # Basis of active flavours
    nbasis = nonzero
    basis = fk.get_flmap()

    # Read the xgrid and the fktable to numpy arrays
    xgrid_flat = fk.get_xgrid()
    fktable_flat = fk.get_sigma()

    # target shape
    if is_hadronic:
        nxsqrt = int(np.sqrt(nx))
        shape_out = (ndata, nbasis, nxsqrt, nxsqrt)
        xgrid = xgrid_flat.reshape(1, nxsqrt)
    else:
        shape_out = (ndata, nbasis, nx)
        xgrid = xgrid_flat.reshape(1, nx)

    # remove padding from the fktable (if any)
    l = len(fktable_flat)
    # get the padding
    pad_position = fk.GetDSz()
    pad_size = fk.GetDSz() - nx * nonzero
    # remove the padding
    if pad_size > 0:
        mask = np.ones(l, dtype=bool)
        for i in range(1, int(l / pad_position) + 1):
            marker = pad_position * i
            mask[marker - pad_size : marker] = False
        fktable_array = fktable_flat[mask]
    else:
        fktable_array = fktable_flat
    # reshape
    fktable = fktable_array.reshape(shape_out)

    dict_out = {
        "ndata": ndata,
        "nbasis": nbasis,
        "nonzero": nbasis,
        "basis": basis,
        "nx": nx,
        "xgrid": xgrid,
        "fktable": fktable,
    }
    return dict_out


def common_data_reader_dataset(dataset_spec):
    """
    Import fktable, common data and experimental data for the given data_name

    Parameters
    ----------
        dataset: validphys representation of a dataset object, a ``validphys.core.DataSetSpec``

    # Arguments:
        - `dataset_c`: c representation of the dataset object
        - `dataset_spec`: python representation of the dataset object

    #Returns:
        - `[dataset_dict]`: a (len 1 list of) dictionary with:
            - dataset: full dataset object from which all data has been extracted
            - hadronic: is hadronic?
            - operation: operation to perform on the fktables below
            - fktables : a list of dictionaries with the following items
                - ndata: number of data points
                - nonzero/nbasis: number of PDF entries which are non zero
                - basis: array containing the information of the non-zero PDF entries
                - nx: number of points in the xgrid (1-D), i.e., for hadronic the total number is nx * nx
                - xgrid: grid of x points
                - fktable: 3/4-D array of shape (ndata, nonzero, nx, (nx))

    instead of the dictionary object that model_gen needs
    """
    dict_fktables = [fk_to_np(fk.load_with_cuts(dataset_spec.cuts)) for fk in dataset_spec.fkspecs]

    # Ask the first fktable for this info
    hadronic = dict_fktables[0]["hadronic"]
    ndata = dict_fktables[0]["ndata"]

    dataset_dict = {
        "fktables": dict_fktables,
        "hadronic": hadronic,
        "operation": dataset_spec.op,
        "name": dataset_spec.name,
        "frac": dataset_spec.frac,
        "ndata": ndata,
    }

    return dataset_dict


def positivity_reader(pos_spec):
    """
    Specific reader for positivity sets
    """
    pos_c = pos_spec.load()
    ndata = pos_c.GetNData()

    parsed_set = [fk_parser(pos_c, pos_c.IsHadronic())]

    pos_sets = [
        {
            "fktables": parsed_set,
            "hadronic": pos_c.IsHadronic(),
            "operation": "NULL",
            "name": pos_spec.name,
            "frac": 1.0,
            "ndata": ndata,
            "tr_fktables": [i["fktable"] for i in parsed_set],
        }
    ]

    positivity_factor = pos_spec.maxlambda

    dict_out = {
        "datasets": pos_sets,
        "trmask": np.ones(ndata, dtype=np.bool),
        "name": pos_spec.name,
        "expdata": np.zeros((1, ndata)),
        "ndata": ndata,
        "positivity": True,
        "lambda": positivity_factor,
        "count_chi2": False,
        "integrability": "INTEG" in pos_spec.name,
    }

    return dict_out

#### new functions ####
def fk_to_np(fkobject):
    """
    Reads a validphys fktable (``validphys.coredata.FKTableData``) and extract
    all necessary information in a clean way, ready to be used with the fitting framework
    """
    #TODO: add this as a method to FKTableData that returns a (tensor, luminosity)
    ndata = fkobject.ndata
    xgrid = fkobject.xgrid
    basis = fkobject.sigma.columns.to_numpy()
    nbasis = len(basis)

    # TODO: I'm sure there's a better way to do this with pandas, please
    # modify it if you know how:
    # The data index needs to be dropped (the data is already cut and is unnecesary sparsing)
    # The flavour index (columns) is correct: only flavours that contribute are included in the tensor
    # The xgrid index instead might need to be filled with 0s or part of the xgrid removed
    # because not all combinations of x-data-flavour contribute

    if fkobject.hadronic:
        # Recover the flavour mapping and make it into a matrix of indices
        ret = np.zeros(14*14, dtype=bool)
        ret[basis] = True
        basis = np.array(np.where(ret.reshape(14, 14))).T.reshape(-1)
        nx = len(xgrid)
        # For non-DIS grids the x1-x2 combinations might be non trivial
        # so I rather not just remove the x i don't like

        # TODO: again, I'm sure there's a better way but the df are not even ordered?
        # get data out of the way
        ns = fkobject.sigma.unstack(["data"], 0)
        # Now let's make sure x1 is complete
        ns = ns.unstack("x2", 0).sort_values("x1").reindex(range(nx), fill_value=0.0)
        # For completeness, the let's ensure the same is true for x2
        ns = ns.stack("x2").unstack("x1", 0).sort_values("x2").reindex(range(nx), fill_value=0.0)
        # Now we have (x2, basis, data, x1) and want (data, basis, x1, x2)
        fktable = np.transpose(ns.values.reshape(nx, nbasis, ndata, nx), (2,1,3,0))
    else:
        # For DIS we can unstack and restack to ensure the ordering makes sense, it should be fast
        full_sigma = fkobject.sigma.unstack(level=["x"],fill_value=0).stack()
        # Turns out that some of the 'x' are actually 0 for all combinations for DIS
        # (and for some reason we kept them!)
        index_x = full_sigma.index.to_frame()["x"].unique()
        xgrid = xgrid[index_x]
        nx = len(xgrid)
        fktable = full_sigma.values.reshape(ndata, nx, nbasis).swapaxes(-2, -1)

    return {
            "ndata": ndata,
            "nbasis": nbasis,
            "nonzero": nbasis,
            "basis": basis,
            "nx": nx,
            "xgrid": xgrid.reshape(1,-1),
            "fktable": fktable,
            "hadronic": fkobject.hadronic
            }

def validphys_group_extractor(group_spec):
    """
    Receives a groupping spec from validphys (most likely and experiment)
    and loops over its content extracting and parsing all information required for the fit

    Parameters
    ----------
        group_spec: Any grouping spec from vp

    Returns
    -------
        parsed_observables: a list of (for now dictionaries) containing the information
                              required to fit the given observable 
    """
    return [common_data_reader_dataset(ds) for ds in group_spec.datasets]
