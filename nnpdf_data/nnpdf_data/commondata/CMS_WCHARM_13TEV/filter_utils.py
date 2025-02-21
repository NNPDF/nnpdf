import logging
import os

from sys_uncertainties import SYS_DEFINITIONS, SYS_UNC_BY_BIN
import yaml

from nnpdf_data.filter_utils.utils import prettify_float, symmetrize_errors

current_dir = os.path.dirname(os.path.abspath(__file__))

yaml.add_representer(float, prettify_float)

MW2 = 80.385**2  # W mass squared in GeV^2
CMSLUMI13 = 2.5  # Luminosity uncertainty in percentage

STAT_LABEL = 'stat_uncorr_unc'


class Extractor:

    def __init__(self, metadata_file, observable, mult_factor=1):
        """
        Extracts kinematics, central data, and uncertainties for a given dataset

        Parameters
        ----------
        metadata_file: str
            Path to the metadata file
        observable: str
            Name of the observable for which the data is extracted. The name
            must be listed in the metadata file.
        mult_factor : float, optional
            Multiplication factor to scale the data. For this dataset it is used
            for a scaling from pb to fb, so a factor 1000.
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

        # Load the (only) table used for this dataset
        table_id = self.metadata["tables"][0]
        with open(f"{current_dir}/rawdata/{table_id}.yaml") as tab:
            self.tab_dict = yaml.safe_load(tab)

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
        [data] = self.tab_dict['independent_variables']
        label = self.metadata['kinematic_coverage']  # ['abs_eta', 'm_W2']
        kinematics = []
        for eta_bin in data['values']:
            abs_eta_min = eta_bin['low']
            abs_eta_max = eta_bin['high']
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
        return kinematics

    def _generate_data_and_unc(self):
        """
        Return a list with central data points and a list with the corresponding
        statistical uncertainties. For this dataset, statistical uncertainties
        are always symmetric.

        The table also provides the corresponding (asymmetric) systematic ucertainty for
        data point. However, this uncertainty is not used as it is preferred to adopt the
        full break-down of the systematic uncertainties. See `_generate_sym_sys_unc`
        """
        logging.info(f"Generating central data for CMS_{self.observable}...")

        [data] = self.tab_dict['dependent_variables']

        # Loop over bins
        dat_central = []
        stat_unc = []
        for rap_bin in data['values']:
            dat_central.append(rap_bin['value'] * self.mult_factor)
            symerror_dict, _asymerror_dict = rap_bin['errors']
            stat_unc.append(symerror_dict['symerror'] * self.mult_factor)

        return dat_central, stat_unc

    def _generate_sym_sys_unc(self):
        """
        The function reads the full break-down of the systematic uncertainties
        as given in the paper. Since such a break-down is not provided in the form of
        a table in HEPData, but rather given as a table in the paper, the list of sources of
        systematic uncertainties is read from an external file (`sys_uncertainties.py`)
        that copies the table in the paper.

        Some of the uncertainties are given in the form of asymmetric uncertainties. These
        asymmetric uncertainties are symmetrized using the usual prescription (see `symmetrize_errors`).

        It returns a list containing a dict for each bin in the absolute rapidity. The keys
        in each dictionary are the names of the sources of uncertainties. The values
        are dicts with keys 'shift', containing the shift from the symmetric prescription, and 'sym_error',
        which is the (symmetrized) value of the uncertainty. Note that the shift is zero if the
        original source of uncertainty is already symmetric.

        Note that uncertainties are given in percentage relative to the central data point
        of the corresponding bin. Moreover, also the shift is a relative value to the central
        data point.
        """
        symmetrized_uncs = []
        for bin in SYS_UNC_BY_BIN:
            unc_dict = {}
            for source in bin:
                if 'asyserror' in source.keys():
                    error_high_low = source['asyserror']
                    plus = error_high_low['high']
                    minus = error_high_low['low']
                    data_delta, sym_error = symmetrize_errors(plus, minus)
                    unc_dict[source['label']] = {'shift': data_delta, 'sym_error': sym_error}
                elif 'syserror' in source.keys():
                    unc_dict[source['label']] = {'shift': 0.0, 'sym_error': source['syserror']}
            symmetrized_uncs.append(unc_dict)
        return symmetrized_uncs

    def _build_unc_definitions(self):
        """
        Build the dictionary containing the definitions of the uncertainties to be
        used in the uncertainty data file.

        The definitions of the systematic uncertainties are given in the
        file `sys_uncertainties.py`.
        """
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
            'type': 'CMSLUMI16',
        }

        # Add systematic uncertainty
        unc_definitions = unc_definitions | SYS_DEFINITIONS

        return unc_definitions

    def generate_data(self):
        '''
        The function collects central data, kinematics, and uncertainties and saves them
        into yaml files.

        The function adds the shifts from the symmetrization prescription to the central
        data points before saving them to the yaml file.

        The systematic uncertainties are given as percentages relative the central data point.
        The absolute value of the uncertainty is obtained from the central data point before
        the shifts are applied.
        '''
        # Get central data, kinematics, and sys uncertainties
        central_data, stat_unc = self._generate_data_and_unc()
        kinematics = self._generate_kinematics()
        symmetrized_sys_uncs = self._generate_sym_sys_unc()

        # Uncertainty definitions
        unc_definitions = self._build_unc_definitions()

        # This loop iterates over the bins of the data.For each bin, it
        # 1) computes the sys_artificial uncertainties, consisting of:
        #    - The effect of symmetrized systematic uncertainties (shift and
        #      sym_error).
        #    - The statistical uncertainty from stat_unc array.
        #    - The luminosity uncertainty.
        # 2) Shifts the central data points central_data[data_idx] to account
        #    for the shift due to  the uncertainty symmetrization
        sys_artificial = []  # Initialize vector of artificial uncertainties
        for data_idx, central_value in enumerate(central_data):
            sys_unc_bin = symmetrized_sys_uncs[data_idx]  # Dict of sys sources for the bin
            shift = 0  # Initialize shift from symmetrization

            # Statistical uncertainty
            unc_dict = {STAT_LABEL: stat_unc[data_idx]}
            # Lmi uncertainty, 0.01 is to convert from percentage to relative value
            unc_dict['corr_lumi_unc'] = central_value * CMSLUMI13 * 0.01

            # Add shift from symmetrization
            for key, value in sys_unc_bin.items():
                # 0.01 is to convert from percentage to relative value
                shift += value['shift'] * 0.01
                unc_dict[key] = value['sym_error'] * central_value * 0.01

            # output of this loop to be saved in the YAML file:
            # 1) list containg uncertainties and
            # 2) central values updated to account for the shift due to symmetization
            sys_artificial.append(unc_dict)
            central_data[data_idx] *= 1.0 + shift

        # Save kinematics into file
        logging.info("Dumping kinematics to file...")
        kinematics_yaml = {'bins': kinematics}
        kins_file_name = self.metadata['kinematics']['file']
        with open(current_dir + "/" + kins_file_name, 'w') as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)
        logging.info("Done!")

        # Save central data into file
        logging.info("Dumping kinematics to file...")
        dat_central_yaml = {'data_central': central_data}
        dat_file_name = self.metadata['data_central']
        with open(current_dir + "/" + dat_file_name, 'w') as file:
            yaml.dump(dat_central_yaml, file, sort_keys=False)
        logging.info("Done!")

        # Save unertainties
        logging.info("Dumping kinematics to file...")
        uncertainties_yaml = {'definitions': unc_definitions, 'bins': sys_artificial}
        unc_file_name = self.metadata['data_uncertainties'][0]
        with open(current_dir + "/" + unc_file_name, 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)
        logging.info("Done!")
