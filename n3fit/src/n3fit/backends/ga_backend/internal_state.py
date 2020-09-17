"""
    Modify the internal state of the backend
"""
from n3fit.backends.ga_backend.MetaModel import MetaModel
from n3fit.backends import select_backend


def set_backend(module):
    """Overwrites the necessary modules and imports from the
    backends module.
    Receives a reference to the module to overwrite.
    """
    # First ensure the backend is tensorflow
    select_backend("Tensorflow")

    # Now set the MetaModel and leave all the rest as the Tensorflow Standard
    setattr(module, "MetaModel", MetaModel)
