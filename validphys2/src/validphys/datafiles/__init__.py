import pathlib
from reportengine.compat import yaml

path_vpdata = pathlib.Path(__file__).parent
path_commondata = pathlib.Path(__file__).with_name('commondata')
path_new_commondata = pathlib.Path(__file__).with_name('new_commondata')
path_theorydb = pathlib.Path(__file__).with_name('theory.db')

# VP should not have access to this file, only to the products
_path_legacy_mapping = path_new_commondata / "dataset_names.yml"
legacy_to_new_mapping = yaml.YAML().load(_path_legacy_mapping)
new_to_legacy_mapping = {}

for k, v in legacy_to_new_mapping.items():
    # Sanitize
    if isinstance(v, str):
        v = {"dataset": v}

    new_dsname = v["dataset"]
    # And fill in the reverse
    new_to_legacy_mapping[new_dsname] = {"dataset": k}
    if "variant" in v:
        new_to_legacy_mapping[new_dsname]["mapping"] = v["variant"]
