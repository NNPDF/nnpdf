import logging
import os

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import matlist_to_matrix, prettify_float, symmetrize_errors

yaml.add_representer(float, prettify_float)

SQRTS = 7000
MW2 = 80.377**2
TABLE_TOKEN = 'Table'
LUMI_UNC = 1.8  # %

# List of systematic uncertainties that shuold
# be considered uncorrelated
SYS_CORR_TOKEN = 'sys_corr_'
SYS_UNCORR_TOKEN = 'sys_uncorr_'


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

        if self.observable == 'WM-YL':
            self.observable_latex = 'W- CJET'
        elif self.observable == 'WP-YL':
            self.observable_latex = 'W+ CBARJET'
        else:
            raise Exception(f'{self.observable} is an unknown observable.')

        # Collect diagonal absoulute uncertainties
        self.diag_unc = self.__collect_diag_unc()

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
        for bin in data['values']:
            pT_min = bin['low']
            pT_max = bin['high']
            kin_bin = {
                label[0]: {'min': pT_min, 'mid': (pT_max + pT_min) / 2, 'max': pT_max},
                label[1]: {'min': None, 'mid': MW2, 'max': None},
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

    def get_table(self, table_id):
        return self.__retrieve_table(table_id)

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

    def generate_data_central(self):
        """
        Same as `generate_kinematics`, but for central data points.
        """
        logging.info(f"Generating central data for ATLAS_{self.observable}...")
        dat_central = []
        table = self.metadata['tables'][0]
        tab_dict = self.__retrieve_table(table)

        # Select the chosen combination
        values = self.__select_bins_by_observable(tab_dict)

        data = [dat['value'] * self.mult_factor for dat in values]
        dat_central = np.concatenate([dat_central, data])
        return dat_central

    def __collect_diag_unc(self):
        """Collect the absolute values of the diagonal uncertainties"""
        table_idx = self.metadata['tables'][0]
        tab_dict = self.__retrieve_table(table_idx)

        # Select the chosen combination
        bins = self.__select_bins_by_observable(tab_dict)

        abs_unc_by_bin = []
        for bin in bins:
            bin_err = bin['errors']
            unc_dict = {
                unc['label']: {'type': list(unc.keys())[1], 'error': list(unc.values())[1]}
                for unc in bin_err
            }
            abs_unc_by_bin.append(unc_dict)
        return abs_unc_by_bin

    def get_diag_unc(self):
        if hasattr(self, 'diag_unc'):
            return self.diag_unc
        else:
            self.diag_unc = self.__collect_diag_unc()
            return self.diag_unc

    def get_unc_def(self):
        if hasattr(self, 'unc_def'):
            return self.unc_def
        else:
            self.unc_def = self.__build_unc_definitions()
            return self.unc_def

    def get_central_data(self):
        if hasattr(self, 'central_data'):
            return self.central_data
        else:
            self.central_data = self.generate_data_central()
            return self.central_data

    def __select_bins_by_observable(self, tab_dict):
        """This dataset colelcts differential xsecs for either W+ and W- in the
        same yaml file. This function returns the part of this yaml file relevant
        for the selected boson sign."""
        values = next(
            (
                head['values']
                for head in tab_dict["dependent_variables"]
                if self.observable_latex in head['header']['name']
            ),
            1,
        )
        if values == 1:
            logging.error(
                f" {self.observable} not found in table under the LaTeX name {self.observable_latex}. The available options are:"
            )
            for head in tab_dict["dependent_variables"]:
                print(f"     - {head['header']['name']}")
            raise Exception()
        else:
            return values

    def __build_unc_definitions(self, variant='default'):
        unc_definitions = {}

        # Statistical uncertainties are always the same
        unc_definitions['stat_uncorr'] = {
            'description': f'Statistical uncertainty',
            'treatment': 'ADD',
            'type': 'UNCORR',
        }

        unc_definitions[SYS_UNCORR_TOKEN + '1'] = {
            'description': f'Uncorrelated systematic uncertainty',
            'treatment': 'MULT',
            'type': 'UNCORR',
        }

        for idx in range(113):
            unc_definitions[SYS_CORR_TOKEN + f'{idx+1}'] = {
                'description': f'Correlated systematic uncertainty idx: {idx+1}',
                'treatment': 'MULT',
                'type': 'CORR',
            }

        # Add Lumi uncertanty
        unc_definitions['LumiUncer'] = {
            'description': f'Integrated luminosity uncertainty',
            'treatment': 'MULT',
            'type': 'ATLASLUMI11',
        }

        if variant == 'sys_10':
            raise ValueError(f"{variant} variant not implemented yet")
        return unc_definitions

    def generate_artifical_unc(self, variant):
        # Collect statistical uncertainties
        diag_uncertainties = self.get_diag_unc()
        stat_unc = [unc['stat']['error'] for unc in diag_uncertainties]
        unc_def = self.get_unc_def()
        data_central = self.get_central_data()

        # Collect systematic uncertainties
        sys_table = self.get_table(self.metadata['tables'][1])['dependent_variables']
        art_sys_unc = []
        for idx_bin, unc_bin in enumerate(sys_table):
            # Get systematic uncertainties
            # NOTE
            # The division 100 is needed if systematic sources are
            # given in percentage in the table (if they are given in percentages)
            bin_sys_unc = [
                value['value'] * data_central[idx_bin] * 0.01 for value in unc_bin['values']
            ]

            # Append Lumi unc
            bin_sys_unc.append(LUMI_UNC * data_central[idx_bin] * 0.01)
            bin_sys_unc = np.concatenate(
                [[stat_unc[idx_bin] * self.mult_factor], bin_sys_unc]
            ).tolist()

            if variant == 'sys_10':
                raise ValueError(f"{variant} variant not implemented yet")

            try:
                assert len(bin_sys_unc) == len(list(unc_def.keys()))
            except AssertionError:
                print(f'problem {len(bin_sys_unc)} != {len(list(unc_def.keys()))}')
            art_sys_unc.append({key: value for key, value in zip(unc_def.keys(), bin_sys_unc)})

        return art_sys_unc

    def generate_data(self, variant='default', save_to_yaml=False, path='./'):
        # Get central data, kinematics, and uncertainties
        central_data = self.get_central_data()
        kinematics = self.generate_kinematics()
        sys_artificial = self.generate_artifical_unc(variant)

        # Uncertainty definitions
        unc_definitions = self.__build_unc_definitions(variant=variant)

        if save_to_yaml:
            # Save kinematics into file
            logging.info("Dumping kinematics to file...")
            kinematics_yaml = {'bins': kinematics}
            with open(path + self.metadata['kinematics']['file'], 'w') as kin_out_file:
                yaml.dump(kinematics_yaml, kin_out_file, sort_keys=False)
            logging.info("Done!")

            # Save central data into file
            logging.info("Dumping kinematics to file...")
            dat_central_yaml = {'data_central': central_data.tolist()}
            file_name = (
                self.metadata['data_central']
                if variant == 'default'
                else self.metadata['variants'][variant]['data_central']
            )
            with open(path + file_name, 'w') as dat_out_file:
                yaml.dump(dat_central_yaml, dat_out_file, sort_keys=False)
            logging.info("Done!")

            # Save unertainties
            logging.info("Dumping kinematics to file...")
            uncertainties_yaml = {'definitions': unc_definitions, 'bins': list(sys_artificial)}
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
