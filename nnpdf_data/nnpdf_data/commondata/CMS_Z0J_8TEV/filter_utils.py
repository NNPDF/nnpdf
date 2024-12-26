import functools
import logging
import os

import numpy as np
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

current_dir = os.path.dirname(os.path.abspath(__file__))

MZ2_LOW = 81.0**2  # GeV2
MZ2_HIGH = 101.0**2  # GeV2
MZ2_MID = (MZ2_LOW + MZ2_HIGH) / 2  # GeV2
CMSLUMI12 = 2.6  # %

ABS_RAP_BINS = [
    {'low': 0.0, 'high': 0.4},
    {'low': 0.4, 'high': 0.8},
    {'low': 0.8, 'high': 1.2},
    {'low': 1.2, 'high': 1.6},
    {'low': 1.6, 'high': 2.0},
]

STAT_ART_LABEL = 'art_corr_unc'
TABLE_TOKEN = 'Table'


class Extractor:
    """
    Extracts kinematics, central data, and uncertainties for a given dataset
    """

    def __init__(self, metadata_file, observable, mult_factor=1.0):
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
        with open(metadata_file) as file:
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
                raise ValueError(f"{observable} is not listed in the metadata file.")

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
        with open(f'{current_dir}/rawdata/{TABLE_TOKEN}{table_id}.yaml') as tab:
            tab_dict = yaml.safe_load(tab)
        return tab_dict

    def _generate_kinematics(self):
        """
        Function that generates the kinematics by taking it from the table with
        measured cross sections. The kinematics are generated in the format of a
        list of dictionaries with the following keys: 'pT', 'abs_eta', 'M_Z2'.
        The values of the keys are dictionaries with the keys 'min', 'mid', and
        'max'.
        """
        logging.info(f"Generating kinematics for CMS_{self.observable}...")

        table = self.metadata["tables"][0]
        tab_dict = self._retrieve_table(table)

        data = tab_dict['independent_variables'][0]
        label = self.metadata['kinematic_coverage']
        kinematics = []
        for rap_bin in ABS_RAP_BINS:
            for pT_bin in data['values']:
                pT_min = pT_bin['low']
                pT_max = pT_bin['high']
                abs_eta_low = rap_bin['low']
                abs_eta_high = rap_bin['high']
                kin_bin = {
                    label[0]: {'min': pT_min, 'mid': (pT_max + pT_min) / 2, 'max': pT_max},
                    label[1]: {
                        'min': abs_eta_low,
                        'mid': (abs_eta_low + abs_eta_high) / 2,
                        'max': abs_eta_high,
                    },
                    label[2]: {'min': MZ2_LOW, 'mid': MZ2_MID, 'max': MZ2_HIGH},
                }
                kinematics.append(kin_bin)

        # Check number of data agrees with metadata
        ndata = len(kinematics)
        if not self.metadata['ndata'] == ndata:
            raise ValueError(
                f"Mismatch in 'ndata': expected {self.metadata['ndata']}, but got {ndata}"
            )

        return kinematics

    def _generate_data_and_unc(self):
        """
        Returns a list with central data points and a list with corresponding uncertainties.
        """
        logging.info(f"Generating central data for CMS_{self.observable}...")
        table = self.metadata['tables'][0]
        tab_dict = self._retrieve_table(table)
        tab_dict = tab_dict['dependent_variables']

        # Loop over kinematic bins
        dat_central = []
        dat_unc = []
        for rap_bin in tab_dict:
            for pt_bin in rap_bin['values']:
                dat_central.append(pt_bin['value'] * self.mult_factor)
                dat_unc.append(pt_bin['errors'][0]['symerror'] * self.mult_factor)

        return dat_central, dat_unc

    def _build_covmat(self):
        '''
        Construct the covarianc matrix from the list of entries provided in HepData.
        '''
        ndata = self.metadata['ndata']
        table_id = self.metadata['tables'][1]
        raw_dict = self._retrieve_table(table_id)

        matlist = [val['value'] for val in raw_dict['dependent_variables'][0]['values']]
        covmat = np.array(matlist).reshape(ndata, ndata)

        if not np.allclose(covmat, covmat.T):
            raise ValueError('Covariance matrix is not symmetric.')

        return covmat

    def _build_unc_definitions(self, variant):
        '''
        Build the dictionary containing the definitions of the uncertainties to
        be used in the uncertainty data file.

        Parameters
        ----------
        variant: str
            Name of the variant to be implemented.

        Return
        ------
        Dict of dicts containing the specifications of each of the
        uncertainties. Each sub-dictionary contains the name of the uncertainty,
        its description, the type, and the treatment. The format is the one used
        in the commondata.
        '''
        unc_definitions = {}

        # Statistical uncertainties are always the same
        for idx in range(self.metadata['ndata']):
            unc_definitions[STAT_ART_LABEL + f'_{idx + 1}'] = {
                'description': f'Artificial uncertainty {idx + 1}, corresponding to a covmat in eigenvector basis',
                'treatment': 'ADD',
                'type': 'CORR',
            }

        # Add lumi uncertainty
        unc_definitions['corr_lumi_unc'] = {
            'description': 'Luminosity uncertainty 2.6%',
            'treatment': 'MULT',
            'type': 'CMSLUMI12',
        }

        if variant == 'sys_10':
            unc_definitions['uncorr_mc_unc'] = {
                'description': 'MC uncertainty',
                'treatment': 'MULT',
                'type': 'UNCORR',
            }

        return unc_definitions

    def generate_data(self, variant='default'):
        '''
        Collect central data, kinematics, and uncertainties and combine them
        in the format used in the commondata.

        Parameters
        ---------
        variant: str
            Name of the dataset variant to generate.
        '''

        # Check if the variant is one of two supported options
        if variant not in ['default', 'sys_10']:
            raise ValueError(f'The variant {variant} is not implemented.')

        # Get central data and kinematics
        central_data, _ = self._generate_data_and_unc()
        kinematics = self._generate_kinematics()

        # Uncertainty definitions
        unc_definitions = self._build_unc_definitions(variant=variant)

        # Get statistical uncertainties. They are represented in a covmat in
        # eigenvector basis, hence they are called "artificial uncertainties".
        # The original covmat can be reconstruted as covat = art_stat.T @ art_stat
        covmat = self._build_covmat()
        eigvals, eigvecs = np.linalg.eig(covmat)
        art_stat = np.sqrt(eigvals) * eigvecs * self.mult_factor

        unc_vals = []  # Initialize vector of uncertainties
        for data_idx, data in enumerate(central_data):
            unc_dict = {}
            for unc_idx, unc_type in enumerate(unc_definitions.keys()):
                if STAT_ART_LABEL in unc_type:
                    # Add statistical uncertainties
                    unc_dict[unc_type] = float(art_stat[data_idx, unc_idx])
                elif unc_type == 'corr_lumi_unc':
                    unc_dict[unc_type] = data * CMSLUMI12 * 0.01
                elif unc_type == 'uncorr_mc_unc':
                    unc_dict[unc_type] = data * 0.01
                else:
                    raise ValueError(f'Uncertainty type {unc_type} is not known.')
            unc_vals.append(unc_dict)

        # Save kinematics into file
        logging.info("Dumping kinematics to file...")
        kinematics_yaml = {'bins': kinematics}
        kins_file_name = self.metadata['kinematics']['file']
        with open(current_dir + '/' + kins_file_name, 'w') as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)
        logging.info("Done!")

        # Save central data into file
        logging.info("Dumping central data to file...")
        dat_central_yaml = {'data_central': central_data}
        data_file_name = self.metadata['data_central']
        with open(current_dir + '/' + data_file_name, 'w') as file:
            yaml.dump(dat_central_yaml, file, sort_keys=False)
        logging.info("Done!")

        # Save unertainties
        logging.info("Dumping uncertainties to file...")
        uncertainties_yaml = {'definitions': unc_definitions, 'bins': unc_vals}
        unc_file_name = (
            self.metadata['data_uncertainties'][0]
            if variant == 'default'
            else self.metadata['variants'][variant]['data_uncertainties'][0]
        )
        with open(current_dir + '/' + unc_file_name, 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)
        logging.info("Done!")
