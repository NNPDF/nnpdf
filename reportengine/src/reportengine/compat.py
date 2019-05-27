"""
compat.py

Fix library compatibility issues.
"""

#To learn about the mess with yaml see this:
#https://bihttps://bitbucket.org/ruamel/yaml/issues/28/consider-making-ruamelyaml-available-as-atbucket.org/ruamel/yaml/issues/28/consider-making-ruamelyaml-available-as-a
try:
    import ruamel_yaml as yaml
except ImportError:
    import ruamel.yaml as yaml

