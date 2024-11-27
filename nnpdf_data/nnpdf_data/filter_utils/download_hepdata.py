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

from hepdata_cli.api import Client, getFilename_fromCd
from hepdata_cli.resilient_requests import resilient_requests
from yaml import safe_load

CLIENT = Client(verbose=True)


class HepDataTables:
    """A commondata class to download the raw HepData tables.

    Parameters:
    -----------
    metadata: str
        a path the metadata of the current data
    """

    def __init__(self, metadata: str) -> None:
        self.metadata = safe_load(pathlib.Path(metadata).read_text())

    def get_hepdata_url(self) -> str:
        # TODO: possibly Check the metadata vs HepData versions here
        hepdata_id = self.metadata["hepdata"]["url"].split("/")[-1]
        # NOTE: `id_list` can support many IDs allowing for downloading
        # mutiple datasets at the same time.
        urls = CLIENT._build_urls(
            id_list=hepdata_id, file_format='yaml', ids='hepdata', table_name=''
        )
        return urls[0]

    def download(self) -> None:
        url = self.get_hepdata_url()
        print(f"Downloading tables from: {url}")
        response = resilient_requests('get', url, allow_redirects=True)

        filename = getFilename_fromCd(response.headers.get('content-disposition'))
        filepath = pathlib.Path().parent
        print(filepath)
        # filepath.mkdir(exist_ok=True)
        raw_path = f"{filepath}/{filename}"

        open(raw_path, 'wb').write(response.content)
        with tarfile.open(raw_path, "r:gz" if raw_path.endswith("tar.gz") else "r:") as tar:
            tar.extractall(path=pathlib.Path(filepath))
        os.remove(raw_path)

        return
