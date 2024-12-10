import logging
import os

import numpy as np
import pandas as pd
from sys_uncertainties import SYS_DEFINITIONS, SYS_UNC_BY_BIN
import yaml

from nnpdf_data.filter_utils.utils import prettify_float, symmetrize_errors

yaml.add_representer(float, prettify_float)

SQRTS = 8000
MW2 = 80.385**2
CMSLUMI13 = 2.5  # %

# List of systematic uncertainties that shuold
# be considered uncorrelated
UNCORR_SYS_UNC = ['UnfoldMCstat', 'UnfoldOtherGen', 'UnfoldReweight']
ART_LABEL = 'art_corr_unc'
STAT_LABEL = 'stat_uncorr_unc'
TABLE = ''

# From Table 1 of the paper
SYS_UNC_by_bin = [{}]


class Extractor:
    """
    Extracts kinematics, central data, and uncertainties for a given dataset

    Parameters
    ----------
    metadata_file: str
      Path to the metadata file
    observable: str
      The name of the observable for which the data is extracted. The name must
      be listed in the metadata file.
    """

    def __init__(self, metadata_file, observable, mult_factor=1):

        # Open metadata and select process
        with open(metadata_file, 'r') as file:
            metadata = yaml.safe_load(file)
            self.metadata = next(
                (
                    md
                    for md in metadata["implemented_observables"]
                    if md['observable_name'] == observable
                ),
                None,
            )
            if self.metadata is None:
                raise Exception(f"{observable} is not listed in the metadata file.")

        # Initialise dict of tables
        self.tables = {}
        self.observable = observable
        self.mult_factor = mult_factor
        self.kin_labels = self.metadata['kinematic_coverage']
        self.ndata = self.metadata['ndata']

        # Collect diagonal absoulute uncertainties
        # self.diag_unc = self.__collect_diag_unc()
        # self.unc_labels = list(self.diag_unc[0].keys())
        # self.unc_labels.pop(0)

    def __retrieve_table(self, table_id):
        """
        Implementation of the lazy loading for the tables. If the table
        is loaded for the first time, it is stored into an internal
        container of the class, so that it will not be loaded each time.

        When called, this functions checks if the table has already been stored
        and, if that is the case, returns the stored table.

        Parameters
        ----------
        table_id: int
          Index that specifies the table

        Return
        ------
        The table specified by `table_id`. If not previously loaded, it is also
        stored into the internal container for future use.
        """
        try:
            table = self.tables[str(table_id)]
        except KeyError:
            logging.debug(
                f'Table {table_id} has not already been used or stored.' f' Storing the table...'
            )
            with open(f'./rawdata/{TABLE}{table_id}.yaml', 'r') as tab:
                tab_dict = yaml.safe_load(tab)
                self.tables[str(table_id)] = tab_dict
                table = tab_dict
        return table

    def __extract_kinematics(self, table: dict):
        """
        Extracts the kinematic variables of the single differential
        distribution given a table.

        For each bin, it computes the max, min, and mid value of the transverse
        momentum of the boson.

        Parameters
        ----------
        table: dict
          Dictionary containing the bins in the transverse momentum
        tab_number: int
          Index to select the range of the second kinematic variable

        Return
        ------
        List of bins containing min, max, and mid values for each of the kinematic
        observables listed in the `kinematic_coverage` of the metadata file.

        """
        data = table['independent_variables'][0]
        label = self.kin_labels
        kinematics = []
        for bin in data['values']:
            abs_eta_min = bin['low']
            abs_eta_max = bin['high']
            kin_bin = {
                label[0]: {
                    'min': abs_eta_min,
                    'mid': (abs_eta_max + abs_eta_min) / 2,
                    'max': abs_eta_max,
                },
                label[1]: {'min': None, 'mid': MW2, 'max': None},
            }
            kinematics.append(kin_bin)
        return kinematics

    def generate_kinematics(self):
        """
        Function that generates the kinematics by looping over all the
        tables specified in the metadata file. The resulting kinematics
        is then saved to a yaml file. It relies on the method
        `__extract_kinematics`.
        """

        logging.info(f"Generating kinematics for ATLAS_{self.observable}...")

        # Initialise kinematics list
        kinematics = []
        ndata = 0
        table = self.metadata["tables"][0]
        tab_dict = self.__retrieve_table(table)
        kin = self.__extract_kinematics(tab_dict)
        kinematics = np.concatenate([kinematics, kin])
        ndata += len(kin)

        # Check number of data agrees with metadata
        try:
            assert self.metadata['ndata'] is not None
            assert self.metadata['ndata'] == ndata
        except AssertionError as e:
            logging.warning(
                f"The number of data in the metafile is either wrong or unspecified."
                f" The correct number is {ndata}. Please, update the metafile."
            )
            return
        return kinematics.tolist()

    def generate_data_and_unc(self, mult_factor=1.0):
        """
        Same as `generate_kinematics`, but for central data points.
        """
        logging.info(f"Generating central data for CMS_{self.observable}...")
        dat_central = []
        stat_unc = []
        asy_sys_unc = []
        table = self.metadata['tables'][0]
        tab_dict = self.__retrieve_table(table)
        tab_dict = tab_dict['dependent_variables'][0]['values']

        # Loop over bins
        for rap_bin in tab_dict:
            dat_central.append(rap_bin['value'] * mult_factor)
            stat_unc.append(rap_bin['errors'][0]['symerror'] * mult_factor)
            asy_sys_unc.append(
                {
                    key: value * mult_factor
                    for key, value in rap_bin['errors'][1]['asymerror'].items()
                }
            )
        return dat_central, stat_unc, asy_sys_unc

    def symmetrized_sys_unc(self):
        """Symmetrise systematic uncertainties. Returns the symmetrized uncertainty
        and the shift to the central data
        """
        symmetrized_uncs = []
        for bin in SYS_UNC_BY_BIN:
            unc_dict = {}
            for source in bin:
                if 'asyserror' in source.keys():
                    error = source['asyserror']
                    plus = error['high']
                    minus = error['low']
                    data_delta, sym_error = symmetrize_errors(plus, minus)
                    unc_dict[source['label']] = {'shift': data_delta, 'sym_error': sym_error}
                elif 'syserror' in source.keys():
                    unc_dict[source['label']] = {'shift': 0.0, 'sym_error': source['syserror']}
            symmetrized_uncs.append(unc_dict)
        return symmetrized_uncs

    def __build_unc_definitions(self, variant='default'):
        unc_definitions = {}

        # Statistical uncertainty
        unc_definitions[STAT_LABEL] = {
            'description': f'Statistical uncertainty',
            'treatment': 'ADD',
            'type': 'UNCORR',
        }

        # Add lumi uncertainty
        unc_definitions['corr_lumi_unc'] = {
            'description': f'Luminosity uncertainty 2.5%',
            'treatment': 'MULT',
            'type': 'CMSLUMI13',
        }

        # Add systematic uncertainty
        unc_definitions = unc_definitions | SYS_DEFINITIONS

        if variant != 'default':
            raise ValueError(f'The variant {variant} is not implemented yet.')

        return unc_definitions

    def generate_data(self, variant='default', save_to_yaml=False, path='./'):
        # Get central data and kinematics
        central_data, stat_unc, _ = self.generate_data_and_unc(self.mult_factor)
        kinematics = self.generate_kinematics()

        # Uncertainty definitions
        unc_definitions = self.__build_unc_definitions(variant=variant)

        sys_artificial = []  # Initialize vector of artificial uncertainties

        symmetrized_sys_uncs = self.symmetrized_sys_unc()
        for data_idx, data in enumerate(central_data):
            shift = 0
            sys_unc_bin = symmetrized_sys_uncs[data_idx]

            # Add shift from symmetrization
            tmp = {}
            for key, value in sys_unc_bin.items():
                shift += value['shift']
                tmp[key] = value['sym_error']

            # Shift central data
            central_data[data_idx] = central_data[data_idx] + shift

            # Statistical uncertainty
            unc_dict = {STAT_LABEL: stat_unc[data_idx]}

            # Lumi uncertainty
            unc_dict['corr_lumi_unc'] = central_data[data_idx] * CMSLUMI13 * 0.01

            # Add systematic uncertainties
            unc_dict = unc_dict | tmp

            sys_artificial.append(unc_dict)

        if save_to_yaml:
            # Save kinematics into file
            logging.info("Dumping kinematics to file...")
            kinematics_yaml = {'bins': kinematics}
            with open(path + self.metadata['kinematics']['file'], 'w') as kin_out_file:
                yaml.dump(kinematics_yaml, kin_out_file, sort_keys=False)
            logging.info("Done!")

            # Save central data into file
            logging.info("Dumping kinematics to file...")
            dat_central_yaml = {'data_central': central_data}
            file_name = self.metadata['data_central']
            with open(path + file_name, 'w') as dat_out_file:
                yaml.dump(dat_central_yaml, dat_out_file, sort_keys=False)
            logging.info("Done!")

            # Save unertainties
            logging.info("Dumping kinematics to file...")
            uncertainties_yaml = {'definitions': unc_definitions, 'bins': sys_artificial}
            file_name = (
                self.metadata['data_uncertainties'][0]
                if variant == 'default'
                else self.metadata['variants'][variant]['data_uncertainties'][0]
            )
            with open(path + file_name, 'w') as dat_out_file:
                yaml.dump(uncertainties_yaml, dat_out_file, sort_keys=False)
            logging.info("Done!")
            return kinematics, central_data, sys_artificial
        else:
            return kinematics, central_data, sys_artificial

    def get_table(self, table_id):
        return self.__retrieve_table(table_id)
