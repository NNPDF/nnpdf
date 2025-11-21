import pathlib

from ._version import __version__
from .commondataparser import parse_new_metadata
from .utils import get_nnpdf_profile
from .validphys_compatibility import legacy_to_new_map

# path_commondata = path_vpdata / "commondata"
THEORY_CARDS_PATH = pathlib.Path(__file__).parent / "theory_cards"


def load_dataset_metadata(dataset_name, variant=None):
    """Given a dataset name, return the metadata"""
    # Compatibility with old nnpdf names, these two lines
    # might disappear at any given point
    if variant is None:
        dataset_name, variant = legacy_to_new_map(dataset_name)

    setname, observable = dataset_name.rsplit("_", 1)

    profile_data_paths = get_nnpdf_profile()["data_path"]

    for data_path in profile_data_paths:
        metadata_file = data_path / setname / "metadata.yaml"
        if metadata_file.exists():
            break
    else:
        raise FileNotFoundError(f"Metadata file for {dataset_name} not found")

    return parse_new_metadata(metadata_file, observable, variant=variant)
