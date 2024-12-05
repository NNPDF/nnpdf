import logging
import os

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import matlist_to_matrix, prettify_float, symmetrize_errors

yaml.add_representer(float, prettify_float)

SQRTS = 8000
MW2 = 80.377**2

# List of systematic uncertainties that shuold
# be considered uncorrelated
UNCORR_SYS_UNC = ['UnfoldMCstat', 'UnfoldOtherGen', 'UnfoldReweight']
STAT_ART_LABEL = 'art_stat_corr'


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

        if self.observable == 'WM-PT':
            self.observable_latex = 'W^-'
        elif self.observable == 'WP-PT':
            self.observable_latex = 'W^+'
        else:
            raise Exception(f'{self.observable} is an unknown observable.')

        # Collect diagonal absoulute uncertainties
        self.diag_unc = self.__collect_diag_unc()
        self.unc_labels = list(self.diag_unc[0].keys())
        self.unc_labels.pop(0)

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
            with open(f'./rawdata/data{table_id}.yaml', 'r') as tab:
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

    def __build_stat_corrmat(self):
        ndata = self.metadata['ndata']
        table_id = self.metadata['tables'][1]
        with open(f'./rawdata/data{table_id}.yaml', 'r') as tab:
            stat = yaml.load(tab, yaml.Loader)
        matlist = [val['value'] for val in stat['dependent_variables'][0]['values']]
        stat_cormat = np.zeros((ndata, ndata))
        for i in range(ndata):
            for j in range(ndata):
                stat_cormat[i, j] = matlist[i + ndata * j]  # Col-major
        return stat_cormat

    def __build_abs_stat_covmat(self):
        """Builds the statistical covmat given the diagonal statistical uncertainties."""
        stat_cormat = self.__build_stat_corrmat()
        # Retrieve statistical diagonal entries
        abs_stat_diag = [float(d['Stat']['error']) for d in self.diag_unc]
        statmat_abs = np.zeros_like(stat_cormat)
        for i in range(stat_cormat.shape[0]):
            for j in range(stat_cormat.shape[1]):
                statmat_abs[i, j] = stat_cormat[i, j] * abs_stat_diag[i] * abs_stat_diag[j]
        return statmat_abs

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

    def symmetrized_sys_unc(self):
        """Symmetrise systematic uncertainties. Returns the symmetrized uncertainty
        and the shift to the central data
        """
        symmetrized_uncs = []
        for bin in self.diag_unc:
            unc_dict = {}
            for source in self.unc_labels:
                if bin[source]['type'] == 'asymerror':
                    error = bin[source]['error']
                    plus = error['plus']
                    minus = error['minus']
                    data_delta, sym_error = symmetrize_errors(plus, minus)
                    unc_dict[source] = {'shift': data_delta, 'sym_error': sym_error}
                elif bin[source]['type'] == 'symerror':
                    unc_dict[source] = {'shift': 0.0, 'sym_error': bin[source]['error']}
            symmetrized_uncs.append(unc_dict)
        return symmetrized_uncs

    # TODO
    # I was not able to reproduce the legacy implementation with this function
    def __cms_sys_unc(self):
        """In the legacy implementation, all systematic uncertainties (even for symmetric uncertanties),
        are treated using the prescription outlined in eq. (6) of https://arxiv.org/pdf/1703.01630.
        Specifically, positive and negative bounds are treated as follows

         C_{ij}^{sys} = \frac{1}{2} \sum_{k'}( C_{jk'}^{+} C_{ik'}^{+} + C_{jk'}^{-} C_{ik'}^{-})

         With this prescription, asymmetric uncertainties are not symmetrized.
        """
        modified_uncertainties = []
        deltas_from_lumi = []
        for bin in self.diag_unc:
            unc_dict = {}
            for source in self.unc_labels:
                if bin[source]['type'] == 'asymerror':
                    error = bin[source]['error']
                    plus = error['plus']
                    minus = error['minus']
                    if source == 'LumiUncert':
                        data_delta, sym_error = symmetrize_errors(plus, minus)
                        unc_dict[f'{source}'] = sym_error
                        data_delta = 0.0
                        deltas_from_lumi.append(data_delta)
                    elif abs(plus) == abs(minus):
                        unc_dict[f'{source}_plus'] = abs(plus) / float(np.sqrt(2))
                        unc_dict[f'{source}_minus'] = 0.0  # minus / float(np.sqrt(2))
                    else:
                        unc_dict[f'{source}_plus'] = plus / float(np.sqrt(2))
                        unc_dict[f'{source}_minus'] = minus / float(np.sqrt(2))
                elif bin[source]['type'] == 'symerror':
                    error = bin[source]['error']
                    if source == 'LumiUncert':
                        unc_dict[f'{source}'] = error
                    elif source in UNCORR_SYS_UNC:
                        unc_dict[f'{source}_plus'] = error / float(np.sqrt(2))
                        unc_dict[f'{source}_minus'] = -error / float(np.sqrt(2))
                    else:
                        unc_dict[f'{source}_plus'] = error / float(np.sqrt(2))
                        unc_dict[f'{source}_minus'] = -error / float(np.sqrt(2))
            modified_uncertainties.append(unc_dict)
        return modified_uncertainties, deltas_from_lumi

    def get_diag_unc(self):
        if hasattr(self, 'diag_unc'):
            return self.diag_unc
        else:
            self.diag_unc = self.__collect_diag_unc()
            return self.diag_unc

    def get_stat_covmat(self):
        return self.__build_abs_stat_covmat()

    def get_stat_cormat(self):
        return self.__build_stat_corrmat()

    def get_abs_stat_covmat(self):
        if hasattr(self, 'abs_stat_covmat'):
            return self.abs_stat_covmat
        else:
            self.abs_stat_covmat = self.__build_abs_stat_covmat()
            return self.abs_stat_covmat

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
                f"{self.observable} not found in table under the LaTeX name {self.observable_latex}. The available options are:"
            )
            for head in tab_dict["dependent_variables"]:
                print(f"     - {head['header']['name']}")
            exit(-1)
        else:
            return values

    def __build_unc_definitions(self, variant='default'):
        unc_definitions = {}

        # Statistical uncertainties are always the same
        for idx in range(self.ndata):
            unc_definitions[STAT_ART_LABEL + f'_{idx + 1}'] = {
                'description': f'Artificial statistical uncertainty {idx + 1}',
                'treatment': 'ADD',
                'type': 'CORR',
            }

        if variant == 'default':
            i = 1
            for label in self.unc_labels:
                if label == 'LumiUncert':
                    unc_type = 'ATLASLUMI12'
                elif label in UNCORR_SYS_UNC:
                    unc_type = 'UNCORR'
                else:
                    unc_type = f'ATLASWJ{i}'
                    i += 1

                unc_definitions[f'{label}'] = {
                    'description': f'Systematic: {label}',
                    'treatment': 'MULT',
                    'type': unc_type,
                }

        elif variant == 'CMS_prescription':
            i = 1
            for label in self.unc_labels:
                if label == 'LumiUncert':
                    unc_definitions[f'{label}'] = {
                        'description': f'Systematic: {label}',
                        'treatment': 'MULT',
                        'type': 'ATLASLUMI12',
                    }
                else:
                    if label in UNCORR_SYS_UNC:
                        unc_definitions[f'{label}_plus'] = {
                            'description': f'Systematic upper unc.: {label}',
                            'treatment': 'MULT',
                            'type': 'UNCORR',
                        }
                        unc_definitions[f'{label}_minus'] = {
                            'description': f'Systematic lower unc.: {label}',
                            'treatment': 'MULT',
                            'type': 'UNCORR',
                        }
                    else:
                        unc_definitions[f'{label}_plus'] = {
                            'description': f'Systematic upper unc.: {label}',
                            'treatment': 'MULT',
                            'type': f'ATLASWJ{i}',
                        }
                        unc_definitions[f'{label}_minus'] = {
                            'description': f'Systematic lower unc.: {label}',
                            'treatment': 'MULT',
                            'type': f'ATLASWJ{i+1}',
                        }
                        i += 2

        elif variant == 'inter_sys':
            raise ValueError(f'Not yet implemented')
        else:
            raise ValueError(f'The variant {variant} is not yet implemented or wrong.')

        return unc_definitions

    def generate_data(self, variant='default', save_to_yaml=False, path='./'):
        # Get central data and kinematics
        central_data = self.generate_data_central()
        kinematics = self.generate_kinematics()

        # Uncertainty definitions
        unc_definitions = self.__build_unc_definitions(variant=variant)

        # Get statistical (artidicial uncertainties)
        stat_covmat = self.__build_abs_stat_covmat()
        eigvals, eigvecs = np.linalg.eig(stat_covmat)
        art_stat = np.sqrt(eigvals) * eigvecs

        sys_artificial = []  # Initialize vector of artificial uncertainties

        if variant == 'default':  # Symmetrized uncs.
            symmetrized_sys_uncs = self.symmetrized_sys_unc()
            shifts = []  # Initialize vector of shifts
            for data_idx, unc in enumerate(symmetrized_sys_uncs):
                shift = 0
                # Add statistical uncertainties
                unc_dict = {
                    STAT_ART_LABEL + f'_{idx + 1}': float(stat_art)
                    for idx, stat_art in enumerate(art_stat[data_idx, :])
                }
                for key, value in unc.items():
                    shift += value['shift']
                    unc_dict[key] = value['sym_error']
                sys_artificial.append(unc_dict)
                shifts.append(shift)
            central_data = central_data + shifts

        elif variant == 'CMS_prescription':
            raise NotImplementedError('This variant has a bug and will not produce any data file.')
            cms_type_uncertainties, deltas = self.__cms_sys_unc()
            central_data = central_data + deltas
            for data_idx, unc in enumerate(cms_type_uncertainties):
                # Add statistical uncertainties
                unc_dict = {
                    STAT_ART_LABEL + f'_{idx + 1}': float(stat_art)
                    for idx, stat_art in enumerate(art_stat[data_idx, :])
                }
                for key, value in unc.items():
                    unc_dict[key] = value
                sys_artificial.append(unc_dict)

        else:
            raise ValueError(f'Variant `{variant}` is not implemented yet.')

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
