"""
    Modify the internal state of the backend
"""
from n3fit.backends.ga_backend.MetaModel import MetaModel


def set_backend(module):
    """Overwrites the necessary modules and imports from the
    backends module.
    Receives a reference to the module to overwrite.
    """
    # Set the MetaModel and leave all the rest as the Tensorflow Standard
    setattr(module, "MetaModel", MetaModel)
