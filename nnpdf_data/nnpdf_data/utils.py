from functools import cache
import os
import pathlib
import pkgutil
import sys

from ruamel.yaml import YAML
from validobj import ValidationError, parse_input

yaml_fast = YAML(typ='safe', pure=False)
yaml_safe = YAML(typ='safe')

# Default path to where the commondata will be found
DEFAULT_PATH_VPDATA = pathlib.Path(__file__).parent / "commondata"
NNPDF_DIR = "NNPDF"


@cache
def get_nnpdf_profile(profile_path=None):
    """Returns the NNPDF profile as a dictionary

    If no ``profile_path`` is provided it will be autodiscovered in the following order:

    1. Environment variable $NNPDF_PROFILE_PATH
    2. ${XDG_CONFIG_HOME}/NNPDF/nnprofile.yaml (usually ~/.config/nnprofile)

    Any value not filled by 1 or 2 will then be filled by the default values
    found within the validphys python package `nnporfile_default.yaml`

    If ``nnpdf_share`` is set to the special key ``RELATIVE_TO_PYTHON``
    the python prefix (``Path(sys.prefix)/"share"/"NNPDF"``) will be used.

    The ``data_path`` key is special in that it will _always_ be filled with
    the location of the nnpdf_data package (albeit with lower priority if there is an entry).
    """

    home_config = pathlib.Path().home() / ".config"
    config_folder = pathlib.Path(os.environ.get("XDG_CONFIG_HOME", home_config)) / NNPDF_DIR

    # Set all default values
    profile_content = pkgutil.get_data("validphys", "nnprofile_default.yaml")
    profile_dict = yaml_safe.load(profile_content)

    # Look at profile path
    if profile_path is None:
        profile_path = os.environ.get("NNPDF_PROFILE_PATH")

    # If profile_path is still none and there is a .config/NNPDF/nnprofile.yaml, read that
    if profile_path is None:
        if (config_nnprofile := config_folder / "nnprofile.yaml").exists():
            profile_path = config_nnprofile
        elif (config_nnprofile := config_folder / "nnprofile.yml").exists():
            profile_path = config_nnprofile

    if profile_path is not None:
        with open(profile_path, encoding="utf-8") as f:
            profile_entries = yaml_safe.load(f)
            if profile_entries is not None:
                profile_dict.update(profile_entries)

    # For the path to the data, loop over the paths we have and check whether they do exist.
    # If any of the paths do not exist, ignore it.
    valid_paths = []
    for dpath_str in profile_dict.get("data_path", []):
        if (dpath := pathlib.Path(dpath_str)).exists():
            valid_paths.append(dpath)
    # And add the default path at the end (i.e., with lower priority)
    valid_paths.append(DEFAULT_PATH_VPDATA)
    profile_dict["data_path"] = valid_paths

    nnpdf_share = profile_dict.get("nnpdf_share")
    if nnpdf_share is None:
        if profile_path is not None:
            raise ValueError(
                f"`nnpdf_share` is not set in {profile_path}, please set it, e.g.: nnpdf_share: `.local/share/NNPDF`"
            )
        raise ValueError(
            "`nnpdf_share` not found in validphys, something is very wrong with the installation"
        )

    if nnpdf_share == "RELATIVE_TO_PYTHON":
        nnpdf_share = pathlib.Path(sys.prefix) / "share" / NNPDF_DIR

    # At this point nnpdf_share needs to be a path to somewhere
    nnpdf_share = pathlib.Path(nnpdf_share)

    # Make sure that we expand any ~ or ~<username>
    nnpdf_share = nnpdf_share.expanduser()

    # Make sure we can either write to this directory or it exists
    try:
        nnpdf_share.mkdir(exist_ok=True, parents=True)
    except PermissionError as e:
        raise FileNotFoundError(
            f"{nnpdf_share} does not exist and you haven't got permissions to create it!"
        ) from e

    # Now read all paths and define them as relative to nnpdf_share (unless given as absolute)
    for var in [
        "results_path",
        "theories_path",
        "validphys_cache_path",
        "hyperscan_path",
        "ekos_path",
        "photons_qed_path",
    ]:
        # if there are any problems setting or getting these variable erroring out is more than justified
        absolute_var = nnpdf_share / pathlib.Path(profile_dict[var]).expanduser()
        profile_dict[var] = absolute_var.absolute().as_posix()

    return profile_dict


def quick_yaml_load(filepath):
    """If libyaml is available, use the C loader to speed up some of the read
    https://pyyaml.org/wiki/LibYAML
    libyaml is available for most linux distributions
    """
    return yaml_fast.load(filepath.read_text(encoding="utf-8"))


def parse_yaml_inp(input_yaml, spec):
    """
    Helper function to parse yaml using the `validobj` library and print
    useful error messages in case of a parsing error.

    https://validobj.readthedocs.io/en/latest/examples.html#yaml-line-numbers
    """
    input_yaml = pathlib.Path(input_yaml)
    inp = quick_yaml_load(input_yaml)
    try:
        return parse_input(inp, spec)
    except ValidationError as e:
        current_exc = e
        # In order to provide a more complete error information, use round trip
        # to read the .yaml file again (insetad of using the CLoader)
        current_inp = YAML(typ="rt").load(input_yaml.open("r", encoding="utf-8"))
        error_text_lines = []
        while current_exc:
            if hasattr(current_exc, 'wrong_field'):
                wrong_field = current_exc.wrong_field
                # Mappings coming from YAML(rt) have an
                # ``lc`` attribute that gives a tuple of
                # ``(line_number, column)`` for a given item in
                # the mapping.
                line = current_inp.lc.item(wrong_field)[0]
                error_text_lines.append(f"Problem processing key at line {line} in {input_yaml}:")
                current_inp = current_inp[wrong_field]
            elif hasattr(current_exc, 'wrong_index'):
                wrong_index = current_exc.wrong_index
                # Similarly lists allow to retrieve the line number for
                # a given item.
                line = current_inp.lc.item(wrong_index)[0]
                current_inp = current_inp[wrong_index]
                error_text_lines.append(
                    f"Problem processing list item at line {line} in {input_yaml}:"
                )
            elif hasattr(current_exc, 'unknown'):
                unknown_lines = []
                for u in current_exc.unknown:
                    unknown_lines.append((current_inp.lc.item(u)[0], u))
                unknown_lines.sort()
                for line, key in unknown_lines:
                    error_text_lines.append(
                        f"Unknown key {key!r} defined at line {line} in {input_yaml}:"
                    )
            error_text_lines.append(str(current_exc))
            current_exc = current_exc.__cause__
        raise ValidationError('\n'.join(error_text_lines)) from e
