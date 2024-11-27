"""
A module to download the HepData tables whenever the commondata
need to be regenerated. As a result, raw data tables are no longer
stored in the repository. It relies on the HepData API for the
downloading of the tables:

https://github.com/HEPData/hepdata-cli/tree/main

TODO: If downloading the raw tables separately for each dataset
turns out to be very slow, a possibility wouldb be to download
all the the raw data tables at once.
"""

import os
import pathlib
import tarfile
from typing import Generator

from hepdata_cli.api import Client, getFilename_fromCd
from hepdata_cli.resilient_requests import resilient_requests
from yaml import safe_load

CLIENT = Client(verbose=True)
COMMONDATA_PATH = pathlib.Path(__file__).parents[1]


class HepDataTables:
    """A commondata class to download the raw HepData tables.

    Parameters
    ----------
    metadata: str
        a path the metadata of the current data
    """

    def __init__(self) -> None:
        self.metadata_files = self._get_metadata_files()

    @staticmethod
    def _get_metadata_files() -> Generator:
        return COMMONDATA_PATH.glob("**/metadata.yaml")

    @staticmethod
    def _get_hepdata_id(metadata_file: pathlib.Path) -> str | None:
        """Get the HepData ID from the metadata for a given dataset.

        Parameters
        ----------
        metadata_file: pathlib.Path
            path to the metadata file

        Returns
        -------
        str | None: returns the HepData ID if the dataset is on HepData
        """
        metadata_yaml = safe_load(pathlib.Path(metadata_file).read_text())
        if "hepdata" in metadata_yaml.keys():
            url = metadata_yaml["hepdata"].get("url", None)
            return url.split("/")[-1] if url is not None else None
        else:
            return None

    @staticmethod
    def _download_data(url: str, dataset_name: str) -> None:
        print(f"Downloading tables from: {url}")
        response = resilient_requests('get', url, allow_redirects=True)

        filename = getFilename_fromCd(response.headers.get('content-disposition'))
        filepath = COMMONDATA_PATH.joinpath(f"commondata/{dataset_name}")
        raw_path = f"{filepath}/{filename}"

        open(raw_path, 'wb').write(response.content)
        with tarfile.open(raw_path, "r:gz" if raw_path.endswith("tar.gz") else "r:") as tar:
            tar.extractall(path=pathlib.Path(filepath))
        os.remove(raw_path)

        return

    def get_hepdata_urls(self) -> tuple[list[str], list[str]]:
        id_list = []
        dataset = []
        for meta in self.metadata_files:
            hep_id = self._get_hepdata_id(metadata_file=meta)
            if hep_id is not None and hep_id.startswith("ins"):
                id_list.append(hep_id)
                dataset.append(str(meta).split("/")[-2])
        hepdata_id = " ".join(id_list)

        urls = CLIENT._build_urls(
            id_list=hepdata_id, file_format='yaml', ids='hepdata', table_name=''
        )
        return urls, dataset

    def download_all(self) -> None:
        urls, dataset_names = self.get_hepdata_urls()
        print("Finished fetching all the URLs.")
        for url, dataname in zip(urls, dataset_names):
            print(f"Downloading: {dataname}")
            self._download_data(url=url, dataset_name=dataname)
        print("Raw data tables for all the datasets have been downloaded.")
        return

    def download_from_server(self) -> None:
        """Download all of the datasets from the NNPDF server."""
        return

    def compress_raw_data(self) -> None:
        """Compress the raw data tables to be uploaded to the NNPDF server."""
        return


def main():
    HepDataTables().download_all()
