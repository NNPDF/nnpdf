try:
    from ._version import __version__
except ModuleNotFoundError:
    raise ModuleNotFoundError("`_version` not found, you might need to reinstall nnpdf (e.g., `pip install .`) if you installed last before nnpdf 4.0.8")
