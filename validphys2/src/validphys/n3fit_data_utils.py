"""
n3fit_data_utils.py

Library of function that read validphys object and make them into clean numpy objects
for its usage for the fitting framework
"""
import numpy as np

def common_data_reader_dataset(dataset_spec):
    """
    Import fktable, common data and experimental data for the given data_name

    Parameters
    ----------
        dataset: validphys representation of a dataset object, a ``validphys.core.DataSetSpec``

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
    parsed_set = [fk_to_np(fk.load_with_cuts(pos_spec.cuts)) for fk in pos_spec.fkspecs]
    hadronic = parsed_set[0]["hadronic"]
    ndata = parsed_set[0]["ndata"]

    pos_sets = [
        {
            "fktables": parsed_set,
            "hadronic": hadronic,
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
    ndata = fkobject.ndata
    xgrid = fkobject.xgrid
    basis = fkobject.luminosity_mapping
    nbasis = len(basis) / 2 if fkobject.hadronic else len(basis)
    nx = len(xgrid)
    fktable = fkobject.get_np_fktable()

    return {
        "ndata": ndata,
        "nbasis": nbasis,
        "nonzero": nbasis,
        "basis": basis,
        "nx": nx,
        "xgrid": xgrid.reshape(1, -1),
        "fktable": fktable,
        "hadronic": fkobject.hadronic,
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
