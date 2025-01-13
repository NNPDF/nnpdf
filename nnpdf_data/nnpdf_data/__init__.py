import pathlib

from ._version import __version__
from .commondataparser import parse_new_metadata
from .validphys_compatibility import legacy_to_new_map, legacy_to_new_mapping, new_to_legacy_map

path_vpdata = pathlib.Path(__file__).parent
path_commondata = path_vpdata / "commondata"
theory_cards = path_vpdata / "theory_cards"


def load_dataset_metadata(dataset_name, variant=None):
    """Given a dataset name, return the metadata"""

    # Compatibility with old nnpdf names, these two lines
    # might disappear at any given point
    if variant is None:
        dataset_name, variant = legacy_to_new_map(dataset_name)

    setname, observable = dataset_name.rsplit("_", 1)
    metadata_file = path_commondata / setname / "metadata.yaml"
    return parse_new_metadata(metadata_file, observable, variant=variant)
