from collections import Counter
import importlib.resources as resources

import pytest

import validphys.scalevariations
from validphys.utils import yaml_safe


def test_unique_theoryid_variations():
    """Check that for each theoryid there is only one set of scale variations in
    scalevariationtheoryids.yaml
    """
    file_path = resources.files(validphys.scalevariations).joinpath("scalevariationtheoryids.yaml")
    with file_path.open("r") as file:
        data = yaml_safe.load(file)
    thids = [k["theoryid"] for k in data["scale_variations_for"]]
    counter = Counter(thids)
    duplicates = [item for item, count in counter.items() if count > 1]
    if duplicates:
        pytest.fail(f"scalevariationtheoryids.yaml multiple entries for theoryIDs: {duplicates}")
