"""
Filter utilities used to build ATLAS Z0J TeV data.

Uncertainty variants:
1) Vanilla: Uncertainties extracted from the source file according the description
provided therein (see below). A luminosity uncertainty is also included.

2) sys_10: On top of 1), this variant includes a 1% uncertainty that takes into
account Monte Carlo uncertainties introduced when grids were first computed.

Source file
------------
The output of each combination is contained within a directory with the same
naming structure as for the input files above. Within these directories are four
files, "*.out , chi2map.dat, tab.dat and xsec1". Again more detailed information
is given in F2AV_README.

The most important file is tab.dat which contains a series of rows, one
corresponding to each bin (after three header rows). The first column is the bin
centre. The next two columns can be ignored. The following columns are then: the
cross section value, the statistical uncertainty (absolute), the systematic
uncertainty uncorrelated between bins (absolute) and the total uncertainty
(absolute) (quadratic sum of statistical + uncorrelated systematics + correlated
systematics) . The following columns detail the orthogonal sources of correlated
systematic uncertainties (relative in percent).
"""

import logging
import os

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

RAP_KINEMATICS_BINS = {
    "Table29": {'range': [0.0, 0.4], 'unc_file': 'ZcombPt_born_m66116_y0004/tab.dat'},
    "Table30": {'range': [0.4, 0.8], 'unc_file': 'ZcombPt_born_m66116_y0408/tab.dat'},
    "Table31": {'range': [0.8, 1.2], 'unc_file': 'ZcombPt_born_m66116_y0812/tab.dat'},
    "Table32": {'range': [1.2, 1.6], 'unc_file': 'ZcombPt_born_m66116_y1216/tab.dat'},
    "Table33": {'range': [1.6, 2.0], 'unc_file': 'ZcombPt_born_m66116_y1620/tab.dat'},
    "Table34": {'range': [2.0, 2.4], 'unc_file': 'ZcombPt_born_m66116_y2024/tab.dat'},
}

M2_KINEMTAICS_BINS = {
    "Table35": {
        'range': np.power([12.0, 20.0], 2).tolist(),
        'unc_file': 'ZcombPt_born_m1220_y0024/tab.dat',
    },
    "Table36": {
        'range': np.power([20.0, 30.0], 2).tolist(),
        'unc_file': 'ZcombPt_born_m2030_y0024/tab.dat',
    },
    "Table37": {
        'range': np.power([30.0, 46.0], 2).tolist(),
        'unc_file': 'ZcombPt_born_m3046_y0024/tab.dat',
    },
    "Table38": {
        'range': np.power([46.0, 66.0], 2).tolist(),
        'unc_file': 'ZcombPt_born_m4666_y0024/tab.dat',
    },
    # "Table39": [66.0, 116.0],
    "Table40": {
        'range': np.power([116.0, 150.0], 2).tolist(),
        'unc_file': 'ZcombPt_born_m116150_y0024/tab.dat',
    },
}

ATLASLUMI12 = 2.8  # percentage
MCUNCERTAINTY = 1.0  # percentage


SQRTS = 8000


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

    def __init__(self, metadata_file, observable, mult_factor):

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

        # Initialise dcit of tables
        self.tables = {}
        self.observable = observable
        self.mult_factor = mult_factor

        # Select the bins for the second kinematic variable
        if observable == 'PT-Y':
            self.kin2_dict = RAP_KINEMATICS_BINS
            self.kin_extra = np.power([66.0, 116.0], 2).tolist()
        elif observable == 'PT-M':
            self.kin2_dict = M2_KINEMTAICS_BINS
            self.kin_extra = [0.0, 2.4]
        else:
            raise Exception(f"Observable {observable} not listed in the metadata file.")
        self.kin_labels = self.metadata['kinematic_coverage']

    def __extract_kinematics(self, table: dict, tab_number: int):
        """
        Extracts the kinematic variables of the double differential
        distribution given a table. Each table represents a bin in the
        second kinematic variable (either rapidity or invariant mass) and
        contains many bins in the transverse momentum of the lepton pair.
        The transverse momentum is taken from `table`, while the second
        kinematic variable is read from the global dicts `RAP_KINEMATICS_BINS`
        or `M_KINEMATICS_BINS` by indicating the `tab_number`.

        For each bin, it computes the max, min, and mid value of the transverse
        momentum and of the second kinematic variable. It also computes the bin for
        a third kinematic variable that specifies an additional kinematic range. For
        instance, in the distribution differential in transverse momentum and rapidity,
        the third kinematic variable is the range in the invariant mass used in the
        measurement. It is understood this range remains the same for all bins. This
        range is selected at initialization, according to the name of the observable.

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
            kin_min = bin['low']
            kin_max = bin['high']
            kin2_bin = self.kin2_dict[f'Table{tab_number}']['range']
            kin_bin = {
                label[0]: {
                    'min': kin2_bin[0],
                    'mid': (kin2_bin[0] + kin2_bin[1]) / 2,
                    'max': kin2_bin[1],
                },
                label[1]: {'min': kin_min, 'mid': (kin_min + kin_max) / 2, 'max': kin_max},
                label[2]: {
                    'min': self.kin_extra[0],
                    'mid': (self.kin_extra[0] + self.kin_extra[1]) / 2,
                    'max': self.kin_extra[1],
                },
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
            with open(f'./rawdata/Table{table_id}.yaml', 'r') as tab:
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

        logging.info(f"Generating kinematics for ATLAS_{self.observable}...")

        # Initialise kinematics list
        kinematics = []
        ndata = 0
        for table in self.metadata["tables"]:
            tab_dict = self.__retrieve_table(table)
            kin = self.__extract_kinematics(tab_dict, table)
            kinematics = np.concatenate([kinematics, kin])
            ndata += len(kin)
        kinematics_yaml = {'bins': kinematics.tolist()}

        # Check number of data agrees with metadata
        try:
            assert self.metadata['ndata'] is not None
            assert self.metadata['ndata'] == ndata
        except AssertionError as e:
            logging.warning(
                f"The number of data in the metafile is either wrong or unspecified."
                f" The correct number is {ndata}. Please, update the metafile."
            )

        # Dump into file
        logging.info("Dumping kinematics to file...")
        with open(self.metadata['kinematics']['file'], 'w') as kin_out_file:
            yaml.dump(kinematics_yaml, kin_out_file, sort_keys=False)
        logging.info("Done!")

    def generate_data_central(self, combination):
        """
        Same as `generate_kinematics`, but for central data points.
        `combination` selects the combination x-sec listed in the table.

        Parameters
        ----------
        combination: str
          Combination x-sec
        """
        logging.info(f"Generating central data for ATLAS_{self.observable}...")
        dat_central = []
        for table in self.metadata['tables']:
            tab_dict = self.__retrieve_table(table)

            # Check if the chosen combination exists
            try:
                assert combination in [
                    head['header']['name'] for head in tab_dict["dependent_variables"]
                ]
            except AssertionError:
                logging.error(f"{combination} is not in table {table}. The available options are:")
                for head in tab_dict["dependent_variables"]:
                    print(f"     - {head['header']['name']}")
                exit(-1)

            # Select the chosen combination
            values = next(
                (
                    head['values']
                    for head in tab_dict["dependent_variables"]
                    if head['header']['name'] == combination
                ),
                None,
            )
            data = [dat['value'] * self.mult_factor for dat in values]
            dat_central = np.concatenate([dat_central, data])

        dat_central_yaml = {'data_central': dat_central.tolist()}

        # Dump into file
        logging.info("Dumping kinematics to file...")
        with open(self.metadata['data_central'], 'w') as dat_out_file:
            yaml.dump(dat_central_yaml, dat_out_file, sort_keys=False)
        logging.info("Done!")

    def generate_uncertainties(self, resource_folder, mc_uncertainity=False):
        """
        Generate the uncertainties and save them into a yaml file. This function
        performs the full break down of the uncertainties as specified in the
        rawdata files. Given that the uncertainties are read from a format that is
        different from the one used to deliver kinematics and central data points,
        this function checks that the ordering of the uncertainties matches the ordering
        of the kinematic variables. The kinematic variables must be previously extracted
        and saved into a yaml file. Failing in that results in an exception.
        """
        logging.info(f"Collecting uncertainties for ATLAS_{self.observable}...")

        MultiIndex = []
        dfs = []
        dirty_flag = False

        for kin in self.kin2_dict.values():
            data_file = kin['unc_file']
            kin_range = kin['range']
            combined_lines = []

            with open(resource_folder + '/' + data_file, 'r') as file:
                lines = file.readlines()

                for i in range(3, len(lines), 2):
                    combined_line = lines[i].strip() + ' ' + lines[i + 1].strip()
                    combined_line = [float(v) for v in combined_line.split()]
                    # Remove useless columns
                    combined_line.pop(1)
                    combined_line.pop(1)
                    # Add lumi. percentage uncertainty
                    combined_line.append(ATLASLUMI12)
                    combined_lines.append(combined_line)

            combined_lines = np.array(combined_lines)
            if not dirty_flag:
                size_of_orth_sys = len(combined_lines[0][5:])
                columns = ['x-sec', 'stat', 'sys_uncorr', 'tot_unc']
                columns = np.concatenate(
                    [
                        columns,
                        [f'sys_corr_{i+1}' for i in range(size_of_orth_sys - 1)],
                        ["sys_lumi_corr"],
                    ]
                )
                if mc_uncertainity:
                    columns = np.append(columns, "sys_mc_uncorr")
                dirty_flag = True

            # Add MC uncertainty
            if mc_uncertainity:
                combined_lines = np.append(
                    combined_lines,
                    [[MCUNCERTAINTY] for _ in range(combined_lines.shape[0])],
                    axis=1,
                )

            dfs.append(pd.DataFrame(combined_lines[:, 1:]))
            MultiIndex.append([(np.mean(kin_range), pt) for pt in combined_lines[:, 0]])

        MultiIndex = np.concatenate(MultiIndex).T
        MultiIndex = pd.MultiIndex.from_arrays(
            MultiIndex, names=[f'{self.kin_labels[0]}_mid', f'{self.kin_labels[1]}_mid']
        )
        df = pd.concat(dfs)
        df.columns = columns
        df = df.set_index(MultiIndex)

        # Multiply relative sys_corr to absolute sys_corr
        for sys_corr in df.columns[4:]:
            df[sys_corr] = df['x-sec'] * df[sys_corr] / 100

        # Multiply all (absolute) values by the multiplicative factor
        df = self.mult_factor * df

        #################################################################################
        # Check that the order of the bins in the unc. dataframe is the same as stored
        # in the kinematics yaml file
        logging.info(f"Checking ordering of bins with {self.metadata['kinematics']['file']}...")
        try:
            assert os.path.exists('./' + self.metadata['kinematics']['file'])
        except AssertionError as e:
            logging.warning(
                "The kinematics file for ATLAS_{self.observable}... has not been "
                "generated yet. Generating it now..."
            )
            self.generate_kinematics()

        with open('./' + self.metadata['kinematics']['file'], 'r') as file:
            kin = yaml.safe_load(file)

            for kin_from_yaml, kin_from_unc in zip(kin['bins'], df.index):
                mid_k1 = (
                    kin_from_yaml[self.kin_labels[0]]['min']
                    + kin_from_yaml[self.kin_labels[0]]['max']
                ) / 2
                mid_pT = (kin_from_yaml['pT']['min'] + kin_from_yaml['pT']['max']) / 2

                try:
                    assert mid_k1 == kin_from_unc[0]
                except AssertionError as e:
                    logging.warning(
                        f"The order of {self.kin_labels[0]} is not the same between the two files. Specifically (HepData) {mid_k1} != {kin_from_unc[0]} (dat)"
                    )

                try:
                    assert mid_pT == kin_from_unc[1]
                except AssertionError as e:
                    logging.warning(
                        f"The order of {self.kin_labels[1]} is not the same between the two files. Specifically (HepData) {mid_pT} != {kin_from_unc[1]} (dat)\n"
                        f"\t  The HepData table is {kin_from_yaml[self.kin_labels[0]]['min']} < {self.kin_labels[0]} < {kin_from_yaml[self.kin_labels[0]]['max']}, "
                        f"and the bin in pT is [{kin_from_yaml['pT']['min']}, {kin_from_yaml['pT']['max']}] with pT_mid = {mid_pT}\n"
                    )

        ####################################

        ###############################
        # Build uncertainty definitions
        unc_definitions = {}
        # Stat unc
        unc_definitions['stat'] = {
            'description': 'Statistical uncorrelated uncertainties',
            'treatment': 'ADD',
            'type': 'UNCORR',
        }

        # Uncorrelated sys unc
        unc_definitions['sys_uncorr'] = {
            'description': 'Systematic uncorrelated uncertainty',
            'treatment': 'MULT',
            'type': 'UNCORR',
        }

        # All the other sources of sys corr unc plus lumi
        for i, unc in enumerate(df.columns[4:]):
            if unc == 'sys_lumi_corr':
                unc_definitions[unc] = {
                    'description': f'Luminosity correlated unc.',
                    'treatment': 'MULT',
                    'type': 'ATLASLUMI12',
                }
            elif unc == 'sys_mc_uncorr':
                unc_definitions[unc] = {
                    'description': f'Monte Carlo uncorrelated uncertainty',
                    'treatment': 'ADD',
                    'type': 'UNCORR',
                }
            else:
                unc_definitions[unc] = {
                    'description': f'Systematic correlated unc. idx: {i+1}',
                    'treatment': 'MULT',
                    'type': 'CORR',
                }

        unc_bins = [
            {key: val for key, val in df.drop(columns=['x-sec', 'tot_unc']).iloc[row].items()}
            for row in range(df.shape[0])
        ]
        unc_yaml = {'definitions': unc_definitions, 'bins': unc_bins}

        # Dump into file
        logging.info("Dumping uncertainties to file...")
        if mc_uncertainity:
            uncertainty_file = self.metadata['variants']['sys_10']['data_uncertainties'][0]
        else:
            uncertainty_file = self.metadata['data_uncertainties'][0]

        with open(uncertainty_file, 'w') as dat_out_file:
            yaml.dump(unc_yaml, dat_out_file, sort_keys=False)
        logging.info("Done!")
