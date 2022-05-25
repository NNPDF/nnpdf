#!/usr/bin/env python

import argparse
import logging
import requests
import sys

from pathlib import Path
from typing import Union

from reportengine.compat import yaml
from reportengine.checks import CheckError


# console = Console()
log = logging.getLogger(__name__)


class DownloadFail(Exception):
    pass


class VersionMismatch(Exception):
    pass


class HepDataConfig:
    
    def __init__(self, metadata_yaml: str, path: Path=None, force: bool=False) -> None:
        """
        Parameters
        ----------
            metadata_yaml: str
                Path to the metadata YAML file
            path: Path
                Path to create the folder in which the tables will be downloaded
            force: bool
                Decide on whether or not overwrite the existing tables
        """

        hep_metadata = self.extract_metadata(metadata_yaml)
        self.url = hep_metadata["url"]
        self.version = int(hep_metadata["version"])
        self.force_download = force
        self.tables = hep_metadata["tables"]

        # Instantiate the dictionary of the scrapped url
        self.hep_webinfo = requests.get(self.url).json()
        self.webversion = self.hep_webinfo["version"]
        self.tables_dict = self.hep_webinfo["hasPart"]
        self.folder = Path().absolute().joinpath("rawdata") if path is None else path

        # Check at initialisation time if local & web versions match. The following is
        # also used to check if the already downloaded HepData tables are outdated. This
        # however assumes any existing tables were downloaded with the version found in
        # the metadata.yaml file (as it should be).
        if self.version != self.webversion and not self.force_download:
            sys.exit(
                f"The required version and the HepData one do not match. The version in the "
                f"metadata is v{self.version} while the one on HepData is v{self.webversion}. "
                f"Please use the flag --force in order to download the new tables."
            )

    def extract_metadata(self, metadata_path: str) -> dict:
        """
        Extract specific information from the metadata, specifically the HepData
        URL and version with the corresponding table number.

        metadata_path: str
            Path to the metadata YAML file
        """
        with open(metadata_path, 'r') as stream_metadata:
            metadata_dic = yaml.safe_load(stream_metadata)
        return metadata_dic["hepdata"]

    def write_tables(self, file: bytes, table_id: Union[int, None]) -> None:
        """
        Write a given table into a file and save it into the disk afterwards.

        Parameters
        ----------
            file: bytes
                HepData table in bytes
            table_id: int
                Denote the `table_id`-th table
        """
        self.folder.mkdir(exist_ok=True)
        ins_id = self.url.split("/")[-1]
        filename = f"HEPData-{ins_id}-v{self.version}-Table_{table_id}"
        # Write the tables as YAML files in rawdata folder
        with open(f"{self.folder}/{filename}.yaml", "wb") as table_yaml:
            table_yaml.write(file)

    def check_downloaded_tables(self, table_numbers: list, vcheck: bool=False) -> None:
        """
        Check if the number of downloaded tables correspond to the list of tables
        in the metadata.yaml file. This function should also be extended to include
        more checks (something much more relevant than this).
        """
        nb_yaml_files_rawdata = sum(1 for _ in self.folder.glob("**/*.yaml")) - 1
        if (nb_yaml_files_rawdata != len(table_numbers)):
            raise DownloadFail("Some of the tables were not downloaded properly.")
        # Check if the version in the hep-metadata and the online one is the same
        if vcheck and Path(f"{self.folder}/hep-metadata-v{self.version}.yaml").is_file():
            with open(f"{self.folder}/hep-metadata-v{self.version}.yaml", "r") as file:
                hep_metadata = yaml.safe_load(file)
            if hep_metadata["version"] != self.webversion and not self.force_download:
                raise VersionMismatch("The local version is different from HepData.")

    def download_heptables(self, table_numbers: list) -> None:
        """
        Download all the HepData tables specified in the metadata.
        """
        # Make sure that the list does not contain duplicates
        table_numbers = list(dict.fromkeys(table_numbers))
        for table in self.tables_dict:
            # Check the current HepData table ID
            table_id = int(table["name"].split()[-1])
            if table_id in table_numbers:
                table_dict = requests.get(table["@id"]).json()
                # Loop over the different table formats. The following surely
                # could be done in a much more efficient way in order to get
                # rid of the recursive loops.
                for url_tab in table_dict["distribution"]:
                    if "YAML" in url_tab["description"]:
                        get_reads_yaml = requests.get(url_tab["contentUrl"])
                        self.write_tables(get_reads_yaml.content, table_id)
                log.info(f"Table {table_id} downloaded and stored properly.")
        # Download the HepData metdata to be used as crosscheck in future.
        # Here, we do not want all the information concerning all the tables
        data = {k:v for k,v in self.hep_webinfo.items() if k != "hasPart"}
        with open(f"{self.folder}/hep-metadata-v{self.version}.yaml", "w") as file:
            yaml.safe_dump(data, file)
        # Crudely check if all the tables have been downloaded successfully.
        # The following check could also be performed much more efficiently.
        self.check_downloaded_tables(table_numbers, vcheck=False)
        print("All the tables have been downloaded and stored properly.")

    def check_hepdata_tables(self) -> None:
        """
        Check if the tables already exist locally. If not, then download all the requiered
        tables. If yes, then perform `check_downloaded_tables`. As mentioned above, this
        should perform more interesting checks (beyond just checking the numbers of tables).
        """
        # Check if the rawdata folder already exists and if so if it contains tables
        if self.folder.exists() and sum(1 for _ in self.folder.glob("**/*.yaml")) > 0:
            self.check_downloaded_tables(self.tables, vcheck=True)
            print("The already downloaded tables match the ones from HepData.")
            return
        self.download_heptables(self.tables)


def argument_parser():
    """Parse input arguments."""
    parser = argparse.ArgumentParser(description="Download & Handle HepData.")
    parser.add_argument("-r", "--runcard", help="Path to YAML metadata run card.")
    parser.add_argument("-f", "--force", action="store_true")
    args = parser.parse_args()

    # TODO: Replace the following with reportengine.checks.make_argcheck
    if not Path(args.runcard).is_file() and not Path(args.runcard).suffix == '.yaml':
        raise CheckError("Invalid runcard. Please use a proper one.")
    return args


def main():
    args = argument_parser()
    hepdata = HepDataConfig(args.runcard, force=args.force)
    hepdata.check_hepdata_tables()


if __name__ == "__main__":
    main()
