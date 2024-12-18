import logging

import numpy as np
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

MW2 = 80.385**2
CMSLUMI13 = 2.5

ART_LABEL = 'art_corr'
STAT_LABEL = 'stat_uncorr'
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

        # Select data with pT > 25 GeV
        tab_dict = tab_dict['dependent_variables'][0]['values']

        # Loop over bins
        for rap_bin in tab_dict:
            dat_central.append(rap_bin['value'] * mult_factor)
            stat_unc.append(rap_bin['errors'][0]['symerror'] * mult_factor)
            asy_sys_unc.append(rap_bin['errors'][1]['symerror'] * mult_factor)
        return dat_central, stat_unc, asy_sys_unc

    def __build_unc_definitions(self):
        unc_definitions = {}

        # Statistical uncertainty
        unc_definitions[STAT_LABEL] = {
            'description': f'Statistical uncertainty',
            'treatment': 'ADD',
            'type': 'UNCORR',
        }

        if self.observable == 'WPWM-RATIO':
            unc_definitions['ART_LABEL'] = {
                'description': f'Correlated systematic uncertainty',
                'treatment': 'MULT',
                'type': 'CORR',
            }
        elif self.observable == 'WPWM-TOT':
            for idx in range(self.ndata):
                unc_definitions[f'{ART_LABEL}_{idx+1}'] = {
                    'description': f'Correlated systematic uncertainty {idx+1}',
                    'treatment': 'ADD',
                    'type': 'CORR',
                }

        return unc_definitions

    def generate_covmat(self, diag_uncs=None):
        table = self.metadata["tables"][1]
        tab_dict = self.__retrieve_table(table)
        matlist = tab_dict['dependent_variables'][0]['values']
        matlist = [d['value'] for d in matlist]
        covmat = np.zeros((self.ndata, self.ndata))
        for i in range(self.ndata):
            for j in range(self.ndata):
                covmat[i, j] = matlist[i + self.ndata * j] * diag_uncs[i] * diag_uncs[j]
        return covmat

    def generate_data(self):
        # Get central data and kinematics
        central_data, stat_unc, sys_unc = self.generate_data_and_unc(self.mult_factor)
        kinematics = self.generate_kinematics()

        # Uncertainty definitions
        unc_definitions = self.__build_unc_definitions()
        sys_artificial = []  # Initialize vector of artificial uncertainties

        if self.observable == 'WPWM-TOT':
            covmat = self.generate_covmat(sys_unc)
            eigvals, eigvecs = np.linalg.eig(covmat)
            art_unc = np.sqrt(eigvals) * eigvecs

            # Loop over bins
            for data_idx, data in enumerate(central_data):
                # Statistical uncertainty
                unc_dict = {STAT_LABEL: stat_unc[data_idx]}
                for sys_idx, art_sys in enumerate(art_unc[data_idx, :]):
                    unc_dict[f'{ART_LABEL}_{sys_idx+1}'] = float(art_sys)
                sys_artificial.append(unc_dict)

        elif self.observable == 'WPWM-RATIO':
            for data_idx, data in enumerate(central_data):
                # Statistical uncertainty
                unc_dict = {STAT_LABEL: stat_unc[data_idx]}

                # Systematic uncertainty
                unc_dict[f'{ART_LABEL}'] = sys_unc[data_idx]
                sys_artificial.append(unc_dict)
        
        # Local path for yaml files
        path='./'

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
        )
        with open(path + file_name, 'w') as dat_out_file:
            yaml.dump(uncertainties_yaml, dat_out_file, sort_keys=False)
        logging.info("Done!")
