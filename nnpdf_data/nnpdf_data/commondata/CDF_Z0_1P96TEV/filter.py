import pathlib

import pandas
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

NB_POINTS = 28
MZ_VALUE = 91.1876  # GeV
SQRT_S = 1_960.0


def load_yaml(table_id: int, version: int = 1) -> dict:
    """Load the HEP data table in yaml format.

    Parameters
    ----------
    table_id: int
        table ID number

    Returns
    -------
    dict:
        ditionary containing the table contents

    """
    filename = f"HEPData-ins856131-v{version}-Table_{table_id}"
    table = pathlib.Path(f"./rawdata/{filename}.yaml")

    return yaml.safe_load(table.read_text())


def get_kinematics(hepdata: dict) -> list:
    """Read the version and list of tables from metadata.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info

    Returns
    -------
    tuple(int, list):
        data version and list of hepdata tables

    """
    rapbins = hepdata["independent_variables"][0]["values"]

    kinematics = []
    for bins in range(NB_POINTS):
        ymin = float(rapbins[bins]["low"])
        ymax = float(rapbins[bins]["high"])
        kin_value = {
            "y": {"min": ymin, "mid": (ymin + ymax) / 2, "max": ymax},
            "m_Z2": {"min": None, "mid": MZ_VALUE**2, "max": None},
            "sqrts": {"min": None, "mid": SQRT_S, "max": None},
        }
        kinematics.append(kin_value)

    return kinematics


def get_data_values(hepdata: dict, indx: int = 0) -> list:
    """Extract the central values from the HepData yaml file.

    Parameters
    ----------
    hepdata: dict
        dictionary containing all data info
    idx: int
        index from which to read the central value, default=0

    Returns
    -------
    list:
        list of dictionaries whose contents are the central values

    """
    central = hepdata["dependent_variables"][indx]["values"]
    return [central[i]["value"] for i in range(NB_POINTS)]


def get_errors() -> dict:
    """Extract the error values from the systematics.dat file.

    Returns
    -------
    pandas.DataFrame:
        dataframe whose contents are the various
        source of uncertainties

    """

    # read the systematics obtained using the c++ script from
    # https://www-cdf.fnal.gov/physics/ewk/2009/dszdy/dszdy_sys.htm
    columns = [
        'y bin',
        'sigma',
        'stat.',
        'lum',
        'B(CC)',
        'B(CP)',
        'B(PP)',
        'CID',
        'PID',
        'CMat',
        'PMat',
        'ZVtx',
        'Trkeff',
        'NoTrk',
        'Tot errors',
    ]

    errors = pd.read_csv("./rawdata/systematics.dat", sep='|', skiprows=3, names=columns)

    return errors


def read_metadata() -> tuple[int, int, list]:
    """Read the version and list of tables from metadata.

    Returns
    -------
    tuple(int, list):
        data version and list of hepdata tables

    """
    metadata = pathlib.Path("./metadata.yaml")
    content = yaml.safe_load(metadata.read_text())

    version = content["hepdata"]["version"]
    nb_datapoints = sum(content["implemented_observables"][0]["npoints"])
    tables = content["implemented_observables"][0]["tables"]

    return version, nb_datapoints, tables


def format_uncertainties(uncs: pandas.DataFrame) -> list:
    """Format the uncertainties to be dumped into the yaml file.

    Parameters
    ----------
    uncs: pandas.DataFrame
        DataFrame containing the various source of uncertainties

    Returns
    -------
    list:
        list of dictionaries whose elements are the various errors

    """

    combined_errors = []
    for _, row in uncs.iterrows():
        error_value = {}
        error_value["stat"] = float(row["stat."])
        for i, sys_i in enumerate(row.iloc[3:-1]):
            error_value[f"sys_corr_{i + 1}"] = sys_i
        combined_errors.append(error_value)

    return combined_errors[:-1]


def dump_commondata(kinematics: list, data: list, errors: pandas.DataFrame) -> None:
    """Function that generates and writes the commondata files.

    Parameters
    ----------
    kinematics: list
        list containing the kinematic values
    data: list
        list containing the central values
    errors: pandas.DataFrame
        DataFrame containing the different errors

    """

    error_definition = {
        "stat": {
            "description": "Uncorrelated statistical uncertainties",
            "treatment": "ADD",
            "type": "UNCORR",
        }
    }

    for i, sys in enumerate(errors.columns[3:-1]):
        error_definition[f"sys_corr_{i + 1}"] = {
            "description": f"Systematic uncertainty {sys}",
            "treatment": "MULT",
            "type": "CORR",
        }

    # update lumi entry
    error_definition['sys_corr_1']['type'] = "CDFLUMI"

    error_definition["stat"] = {
        "description": "Uncorrelated statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    errors_formatted = format_uncertainties(errors)

    with open("data_ZRAP.yaml", "w") as file:
        yaml.dump({"data_central": data}, file, sort_keys=False)

    with open("kinematics_ZRAP.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open("uncertainties_ZRAP.yaml", "w") as file:
        yaml.dump(
            {"definitions": error_definition, "bins": errors_formatted}, file, sort_keys=False
        )


def main_filter() -> None:
    """Main driver of the filter that produces commmondata.

    There are four main different sources of uncertainties.

    1. Statistical uncertainties: ADD, UNCORR

    2. Correlated Systematic uncertainties: MULT, CORR:
        Obtained from a c++ script provided with the experimental paper 0908.3914

    4. Luminosity Systematic uncertainties: MULT, CDFLUMI

    """

    yaml_content = load_yaml(table_id=2, version=1)

    kinematics = get_kinematics(yaml_content)
    data_central = get_data_values(yaml_content)
    uncertainties = get_errors()

    # correlations from https://inspirehep.net/literature/806697
    # compile using
    # g++ -c error_propagator_g++_032610.C
    # g++ error_propagator_g++_032610.o -o systematics

    # Generate all the necessary files
    dump_commondata(kinematics, data_central, uncertainties)

    return


if __name__ == "__main__":
    main_filter()
