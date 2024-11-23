import pathlib

from ruamel.yaml import YAML
from validobj import ValidationError, parse_input
import yaml

yaml_rt = YAML(typ="rt")
try:
    # If libyaml is available, use the C loader
    Loader = yaml.CLoader
except AttributeError:
    # fallback to the slow loader
    Loader = yaml.Loader


def quick_yaml_load(filepath):
    """If libyaml is available, use the C loader to speed up some of the read
    https://pyyaml.org/wiki/LibYAML
    libyaml is available for most linux distributions
    """
    return yaml.load(filepath.read_text(encoding="utf-8"), Loader=Loader)


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
        current_inp = yaml_rt.load(input_yaml.open("r", encoding="utf-8"))
        error_text_lines = []
        while current_exc:
            if hasattr(current_exc, 'wrong_field'):
                wrong_field = current_exc.wrong_field
                # Mappings compping from yaml_rt have an
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
