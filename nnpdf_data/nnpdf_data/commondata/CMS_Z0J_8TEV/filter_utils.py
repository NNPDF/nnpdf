import logging
import os

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import prettify_float, symmetrize_errors

yaml.add_representer(float, prettify_float)

SQRTS = 8000
MZ2_low = 81.0**2  # GeV2
MZ2_high = 101.0**2  # GeV2
MZ2_mid = (MZ2_low + MZ2_high) * 0.5  # GeV2
CMSLUMI12 = 2.6  # %

ABS_RAP_BINS = [
    {'low': 0.0, 'high': 0.4},
    {'low': 0.4, 'high': 0.8},
    {'low': 0.8, 'high': 1.2},
    {'low': 1.2, 'high': 1.6},
    {'low': 1.6, 'high': 2.0},
]

# List of systematic uncertainties that shuold
# be considered uncorrelated
UNCORR_SYS_UNC = ['UnfoldMCstat', 'UnfoldOtherGen', 'UnfoldReweight']
STAT_ART_LABEL = 'art_corr_unc'
TABLE_TOKEN = 'Table'


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

    def __init__(self, metadata_file, observable, mult_factor=1.0):

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

    def __extract_kinematics(self, table: dict, tab_number: int):
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
        for rap_bin in ABS_RAP_BINS:
            for bin in data['values']:
                pT_min = bin['low']
                pT_max = bin['high']
                abs_eta_low = rap_bin['low']
                abs_eta_high = rap_bin['high']
                kin_bin = {
                    label[0]: {'min': pT_min, 'mid': (pT_max + pT_min) / 2, 'max': pT_max},
                    label[1]: {
                        'min': abs_eta_low,
                        'mid': (abs_eta_low + abs_eta_high) * 0.5,
                        'max': abs_eta_high,
                    },
                    label[2]: {'min': MZ2_low, 'mid': MZ2_mid, 'max': MZ2_high},
                }
                kinematics.append(kin_bin)
        return kinematics

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
            with open(f'./rawdata/{TABLE_TOKEN}{table_id}.yaml', 'r') as tab:
                tab_dict = yaml.safe_load(tab)
                self.tables[str(table_id)] = tab_dict
                table = tab_dict
        return table

    def generate_kinematics(self):
        """
        Function that generates the kinematics by looping over all the
        tables specified in the metadata file. The resulting kinematics
        is then saved to a yaml file. It relies on the method
        `__extract_kinematics`.
        """

        logging.info(f"Generating kinematics for CMS_{self.observable}...")

        # Initialise kinematics list
        kinematics = []
        ndata = 0
        table = self.metadata["tables"][0]
        tab_dict = self.__retrieve_table(table)
        kin = self.__extract_kinematics(tab_dict, table)
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
        dat_unc = []
        table = self.metadata['tables'][0]
        tab_dict = self.__retrieve_table(table)
        tab_dict = tab_dict['dependent_variables']

        # Loop over rap bins
        for rap_bin in tab_dict:
            for pt_bin in rap_bin['values']:
                dat_central.append(pt_bin['value'] * mult_factor)
                dat_unc.append(pt_bin['errors'][0]['symerror'] * mult_factor)
        return dat_central, dat_unc

    def build_covmat(self):
        ndata = self.metadata['ndata']
        table_id = self.metadata['tables'][1]
        with open(f'./rawdata/{TABLE_TOKEN}{table_id}.yaml', 'r') as tab:
            raw_dict = yaml.load(tab, yaml.Loader)
        matlist = [val['value'] for val in raw_dict['dependent_variables'][0]['values']]
        covmat = np.zeros((ndata, ndata))
        for i in range(ndata):
            for j in range(ndata):
                covmat[i, j] = matlist[i + ndata * j]  # Col-major
        return covmat

    def __build_unc_definitions(self, variant='default'):
        unc_definitions = {}

        # Statistical uncertainties are always the same
        for idx in range(self.ndata):
            unc_definitions[STAT_ART_LABEL + f'_{idx + 1}'] = {
                'description': f'Artificial uncertainty {idx + 1}',
                'treatment': 'ADD',
                'type': 'CORR',
            }

        # Add lumi uncertainty
        unc_definitions['corr_lumi_unc'] = {
            'description': f'Luminosity uncertainty 2.6%',
            'treatment': 'MULT',
            'type': 'CMSLUMI12',
        }

        if variant == 'sys_10':
            unc_definitions['uncorr_mc_unc'] = {
                'description': f'MC uncertainty',
                'treatment': 'MULT',
                'type': 'UNCORR',
            }
        elif variant != 'default':

            raise ValueError(f'The variant {variant} is not implemented yet.')

        return unc_definitions

    def generate_data(self, variant='default', save_to_yaml=False, path='./'):
        # Get central data and kinematics
        central_data, _ = self.generate_data_and_unc(self.mult_factor)
        kinematics = self.generate_kinematics()

        # Uncertainty definitions
        unc_definitions = self.__build_unc_definitions(variant=variant)

        # Get statistical (artidicial uncertainties)
        covmat = self.build_covmat()
        eigvals, eigvecs = np.linalg.eig(covmat)
        art_stat = np.sqrt(eigvals) * eigvecs * self.mult_factor
        sys_artificial = []  # Initialize vector of artificial uncertainties

        for data_idx, data in enumerate(central_data):
            unc_dict = {}
            # Add artificial uncertainties
            for unc_idx, unc_type in enumerate(unc_definitions.keys()):
                if STAT_ART_LABEL in unc_type:
                    # Add statistical uncertainties
                    unc_dict[unc_type] = float(art_stat[data_idx, unc_idx])
                elif unc_type == 'corr_lumi_unc':
                    unc_dict[unc_type] = data * CMSLUMI12 * 0.01
                elif unc_type == 'uncorr_mc_unc' and variant == 'sys_10':
                    unc_dict[unc_type] = data * 0.01
                else:
                    raise ValueError(f'Uncertainty type {unc_type} is not known.')
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

    # Getters
    def get_table(self, table_id):
        return self.__retrieve_table(table_id)

    def get_diag_unc(self):
        if hasattr(self, 'diag_unc'):
            return self.diag_unc
        else:
            _, self.diag_unc = self.generate_data_and_unc()
            return self.diag_unc

    def get_covmat(self):
        return self.build_covmat()
