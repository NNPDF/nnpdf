import logging
import os

import numpy as np
import pandas as pd
import yaml

from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

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

        if self.observable == 'WM-PT':
            self.observable_latex = 'W^-'
        elif self.observable == 'WP-PT':
            self.observable_latex = 'W^+'
        else:
            raise Exception(f'{self.observable} is an unknown observable.')

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
            kin_bin = {label[0]: {'min': pT_min, 'mid': (pT_max + pT_min) / 2, 'max': pT_max}}
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

    def generate_data_central(self):
        """
        Same as `generate_kinematics`, but for central data points.
        """
        logging.info(f"Generating central data for ATLAS_{self.observable}...")
        dat_central = []
        for table in self.metadata['tables']:
            tab_dict = self.__retrieve_table(table)

            # Select the chosen combination
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
                    f"{self.observable} not found in table {table} under the LaTeX name {self.observable_latex}. The available options are:"
                )
                for head in tab_dict["dependent_variables"]:
                    print(f"     - {head['header']['name']}")
                exit(-1)

            data = [dat['value'] * self.mult_factor for dat in values]
            dat_central = np.concatenate([dat_central, data])

        dat_central_yaml = {'data_central': dat_central.tolist()}

        # Dump into file
        logging.info("Dumping kinematics to file...")
        with open(self.metadata['data_central'], 'w') as dat_out_file:
            yaml.dump(dat_central_yaml, dat_out_file, sort_keys=False)
        logging.info("Done!")
