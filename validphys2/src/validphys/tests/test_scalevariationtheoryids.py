import importlib.resources as resources

from ruamel.yaml import YAML

import validphys.scalevariations

yaml = YAML()


def test_unique_theoryid_variations():
    """Check that for each theoryid there is only one set of scale variations in
    scalevariationtheoryids.yaml
    """
    file_path = resources.files(validphys.scalevariations).joinpath("scalevariationtheoryids.yaml")
    with file_path.open("r") as file:
        data = yaml.load(file)
    thids = [k["theoryid"] for k in data["scale_variations_for"]]
    assert len(thids) == len(set(thids))
