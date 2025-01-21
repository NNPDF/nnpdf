import functools
import logging
import os

import numpy as np
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

MW2 = 80.385**2

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ART_LABEL = 'art_corr'
STAT_LABEL = 'stat_uncorr'
TABLE_TOKEN = 'Table'


class Extractor:

    def __init__(self, metadata_file, observable, mult_factor=1):
        """
        Parameters
        ----------
        metadata_file: str
            Path to the metadata file
        observable: str
            The name of the observable for which the data is extracted. The name
            must be listed in the metadata file.
        mult_factor: float
            Multiplication factor to apply to the central data points. This is
            useful to convert the data in the metadata file to the desired
            units.
        """
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

        self.observable = observable
        self.mult_factor = mult_factor

    @functools.cache
    def _retrieve_table(self, table_id):
        """
        Implementation of the loading for the table.

        Parameters
        ----------
        table_id: int
          Index that specifies the table.

        Return
        ------
        The table specified by `table_id`.
        """
        with open(f'{CURRENT_DIR}/rawdata/{TABLE_TOKEN}{table_id}.yaml') as tab:
            tab_dict = yaml.safe_load(tab)
        return tab_dict

    def _generate_kinematics(self):
        """
        The function generates the kinematics by reading and processing it from
        the referenced table. Kinematics is processed in the format of a list of
        dictionaries. The keys in each dictionaries specify the label (i.e. name)
        for the kinematic variables. For this dataset, they are 'abs_eta' and 'm_W2'.
        The labels are taken from the matadata file. The corresponding values are
        'min', 'mid', and 'max'.

        For this dataset, 'm_W2' is used in the computation of the (x,Q2)-map and
        does not have any active role in the fit. For that reason, every bin has the
        same value. Moreover, only the mid value is used.
        """
        logging.info(f"Generating kinematics for CMS_{self.observable}...")

        table_ID = self.metadata["tables"][0]
        tab_dict = self._retrieve_table(table_ID)

        data = tab_dict['independent_variables'][0]
        label = self.metadata['kinematic_coverage']
        kinematics = []
        for eta_bin in data['values']:
            abs_eta_max = eta_bin['high']
            abs_eta_min = eta_bin['low']
            kin_bin = {
                label[0]: {
                    'min': abs_eta_min,
                    'mid': (abs_eta_max + abs_eta_min) / 2,
                    'max': abs_eta_max,
                },
                label[1]: {'min': None, 'mid': MW2, 'max': None},
            }
            kinematics.append(kin_bin)

        # Check number of data agrees with metadata
        ndata = len(kinematics)
        if not self.metadata['ndata'] == ndata:
            raise ValueError(
                f"Mismatch in 'ndata': expected {self.metadata['ndata']}, but got {ndata}"
            )
        self.ndata = ndata
        return kinematics

    def _generate_data_and_unc(self):
        """
        Return a list with central data points and two additional lists with the corresponding
        statistical and systematic uncertainties. For this dataset, uncertainties are always
        symmetric. Uncertainties are given as absolute values.

        Note that, for the total x-sec, the correlation matrix is provided. The corresponding
        covariance matrix is constructed in `_generate_covmat`.
        """
        logging.info(f"Generating central data for CMS_{self.observable}...")
        dat_central = []
        stat_unc = []
        asy_sys_unc = []
        table_ID = self.metadata['tables'][0]
        tab_dict = self._retrieve_table(table_ID)

        # Select data with pT > 25 GeV
        tab_dict = tab_dict['dependent_variables'][0]['values']

        # Loop over bins
        for rap_bin in tab_dict:
            dat_central.append(rap_bin['value'] * self.mult_factor)
            stat_unc.append(rap_bin['errors'][0]['symerror'] * self.mult_factor)
            asy_sys_unc.append(rap_bin['errors'][1]['symerror'] * self.mult_factor)
        return dat_central, stat_unc, asy_sys_unc

    def _build_unc_definitions(self):
        """
        Build the dictionary containing the definitions of the uncertainties to be
        used in the uncertainty data file.
        """
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

    def _generate_covmat(self, diag_uncs):
        """
        Generate the covariance matrix for the total x-sec. This function requires
        the diagonal systematic uncertainties as argument. The diagonal uncertainties
        are used to construct the covariance matrix from the correlation matrix stored
        in the HepData table.

        Note that such a correlation matrix exists for the total x-sec only, while the
        ratio observable does not provide this information.
        """
        if not self.observable == 'WPWM-TOT':
            raise ValueError(
                "The construction of the covariance matrix is defined for the total x-sec only."
            )
        table_ID = self.metadata["tables"][1]
        tab_dict = self._retrieve_table(table_ID)
        matlist = tab_dict['dependent_variables'][0]['values']
        matlist = [d['value'] for d in matlist]
        covmat = np.zeros((self.ndata, self.ndata))
        for i in range(self.ndata):
            for j in range(self.ndata):
                covmat[i, j] = matlist[i + self.ndata * j] * diag_uncs[i] * diag_uncs[j]
        return covmat

    def generate_data(self):
        """
        The function collects central data, kinematics, and uncertainties ans save them
        into yaml files.

        The systematic uncertainties are given as percentages relative the central data point.
        The absolute value of the uncertainty is obtained from the central data point before
        the shifts are applied.
        """
        # Get central data and kinematics
        central_data, stat_unc, sys_unc = self._generate_data_and_unc()
        kinematics = self._generate_kinematics()

        # Uncertainty definitions
        unc_definitions = self._build_unc_definitions()
        sys_artificial = []  # Initialize vector of artificial uncertainties

        if self.observable == 'WPWM-TOT':
            # Generate covmat and perform eigen decomposition
            covmat = self._generate_covmat(sys_unc)
            eigvals, eigvecs = np.linalg.eig(covmat)
            art_unc = np.sqrt(eigvals) * eigvecs

            # Loop over bins
            for data_idx in range(len(central_data)):
                # Statistical uncertainty
                unc_dict = {STAT_LABEL: stat_unc[data_idx]}

                # Artificial systematic uncertainties
                for sys_idx, art_sys in enumerate(art_unc[data_idx, :]):
                    unc_dict[f'{ART_LABEL}_{sys_idx+1}'] = float(art_sys)

                # Append to list
                sys_artificial.append(unc_dict)

        elif self.observable == 'WPWM-RATIO':
            for data_idx in range(len(central_data)):
                # Statistical uncertainty
                unc_dict = {STAT_LABEL: stat_unc[data_idx]}

                # Systematic uncertainty
                unc_dict[f'{ART_LABEL}'] = sys_unc[data_idx]
                sys_artificial.append(unc_dict)

        # Save kinematics into file
        logging.info("Dumping kinematics to file...")
        kinematics_yaml = {'bins': kinematics}
        kins_file_name = self.metadata['kinematics']['file']
        with open(CURRENT_DIR + '/' + kins_file_name, 'w') as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)
        logging.info("Done!")

        # Save central data into file
        logging.info("Dumping kinematics to file...")
        dat_central_yaml = {'data_central': central_data}
        data_file_name = self.metadata['data_central']
        with open(CURRENT_DIR + '/' + data_file_name, 'w') as file:
            yaml.dump(dat_central_yaml, file, sort_keys=False)
        logging.info("Done!")

        # Save unertainties
        logging.info("Dumping kinematics to file...")
        uncertainties_yaml = {'definitions': unc_definitions, 'bins': sys_artificial}
        unc_file_name = self.metadata['data_uncertainties'][0]
        with open(CURRENT_DIR + '/' + unc_file_name, 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)
        logging.info("Done!")
