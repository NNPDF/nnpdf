import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.legacy_jets_utils import (
    TABLE_TO_RAPIDITY_ATLAS_1JET_7TEV,
    fill_df_ATLAS_1JET_7TEV,
)
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

AVAILABLE_VARIANTS = ['nominal', 'weaker', 'stronger']


def observable_selector(metadata: dict, observable: str):
    """
    Select the metadata for the given observable.
    """
    for obs in metadata["implemented_observables"]:
        if obs["observable_name"] == observable:
            return obs
    raise ValueError(f"Observable {observable} not found in metadata.")


def filter_ATLAS_1JET_7TEV_data_kinematics(observable: str = 'PTY-R06'):
    """
    Write kinematic values in the kinematics.yaml file.
    """
    with open("metadata.yaml") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    obs_metadata = observable_selector(metadata, observable)
    tables = obs_metadata["tables"]
    kinematics_filename = obs_metadata["kinematics"]["file"]

    kin = []
    for table in tables:
        hepdata_table = f"rawdata/HEPData-ins1325553-v{version}_table{table}.yaml"

        with open(hepdata_table) as file:
            input = yaml.safe_load(file)

        # rapidity
        rapidity_interval = TABLE_TO_RAPIDITY_ATLAS_1JET_7TEV[table]
        rap = {}
        rap['min'], rap['max'] = rapidity_interval[0], rapidity_interval[1]
        rap['mid'] = 0.5 * (rap['min'] + rap['max'])

        # center of mass energy
        sqrts = float(input['dependent_variables'][0]['qualifiers'][4]['value'])

        # transverse momentum
        jet_kt_bins = input['independent_variables'][0]['values']
        KT = {}
        for kt in jet_kt_bins:
            KT['min'], KT['max'] = kt['low'], kt['high']
            KT['mid'] = float(f"{0.5 * (kt['low'] + kt['high']):.3f}")

            kin_value = {
                'y': {'min': rap['min'], 'mid': rap['mid'], 'max': rap['max']},
                'pT': {'min': KT['min'], 'mid': KT['mid'], 'max': KT['max']},
                'sqrts': {'min': None, 'mid': sqrts, 'max': None},
            }

            kin.append(kin_value)

    kinematics_yaml = {"bins": kin}

    with open(kinematics_filename, "w") as file:
        yaml.dump(kinematics_yaml, file, sort_keys=False)


def filter_ALTAS_1JET_7TEV_data_uncertainties(observable: str = 'PTY-R06', variant='nominal'):
    """
    Write uncertainties in the uncertainties.yaml file.
    """

    if variant not in AVAILABLE_VARIANTS:
        raise ValueError(
            f"Variant {variant} not recognized. Available variants: {AVAILABLE_VARIANTS}"
        )

    with open("metadata.yaml") as file:
        metadata = yaml.safe_load(file)

    version = metadata["hepdata"]["version"]
    obs_metadata = observable_selector(metadata, observable)
    tables = obs_metadata["tables"]

    # get df of uncertainties
    dfs = []
    cvs = []
    for table in tables:
        # uncertainties dataframe
        df, cv = fill_df_ATLAS_1JET_7TEV(table, version, variant)
        dfs.append(df)
        cvs.append(cv)

    df_unc = pd.concat([df for df in dfs], axis=0)
    cvs = np.concatenate(cvs, axis=0)

    # statistical errors fully uncorrelated
    stat_errors = df_unc["stat"].to_numpy()

    # luminosity errors
    lum_errors = df_unc["sys_lumi"].to_numpy()

    A_corr = df_unc.drop(["stat", "sys_lumi"], axis=1).to_numpy()

    # Error definitions
    error_definition = {
        f"{col}": {
            "description": f"correlated systematic {col}",
            "treatment": "MULT",
            "type": "CORR",
        }
        for col in df_unc.drop(["stat", "sys_lumi"], axis=1).columns
    }

    error_definition["luminosity_uncertainty"] = {
        "description": "luminosity uncertainty",
        "treatment": "MULT",
        "type": "ATLASLUMI14",
    }

    error_definition["statistical_uncertainty"] = {
        "description": "statistical uncertainty",
        "treatment": "MULT",
        "type": "UNCORR",
    }

    # store error in dict
    error = []
    for n in range(A_corr.shape[0]):
        error_value = {}
        for col, m in zip(
            df_unc.drop(["stat", "sys_lumi"], axis=1).columns, range(A_corr.shape[1])
        ):
            error_value[f"{col}"] = float(A_corr[n, m])

        error_value["luminosity_uncertainty"] = float(lum_errors[n])
        error_value["statistical_uncertainty"] = float(stat_errors[n])
        error.append(error_value)

    uncertainties_yaml = {"definitions": error_definition, "bins": error}
    if variant == 'nominal':
        data_filename = obs_metadata['data_central']
        uncertainties_filename = obs_metadata['data_uncertainties'][0]
    else:
        data_filename = obs_metadata['variants'][variant]['data_central']
        uncertainties_filename = obs_metadata['variants'][variant]['data_uncertainties'][0]

    with open(uncertainties_filename, "w") as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)

    data_central_yaml = {"data_central": cvs.tolist()}
    with open(data_filename, "w") as file:
        yaml.dump(data_central_yaml, file, sort_keys=False)


if __name__ == "__main__":
    filter_ATLAS_1JET_7TEV_data_kinematics(observable='PTY-R06')
    filter_ALTAS_1JET_7TEV_data_uncertainties(observable='PTY-R06', variant='nominal')
    filter_ALTAS_1JET_7TEV_data_uncertainties(observable='PTY-R06', variant='weaker')
    filter_ALTAS_1JET_7TEV_data_uncertainties(observable='PTY-R06', variant='stronger')
    filter_ALTAS_1JET_7TEV_data_uncertainties(observable='PTY-R04', variant='nominal')
    filter_ALTAS_1JET_7TEV_data_uncertainties(observable='PTY-R04', variant='weaker')
    filter_ALTAS_1JET_7TEV_data_uncertainties(observable='PTY-R04', variant='stronger')
