# Make the chosen backend available from n3fit.backends
# These includes:
#   set_initial_state, clear_backend_state
#   meta classes: MetaLayer, MetaModel
#   base_layers: Input, Lambda, base_layer_selector, regularizer_selector, Concatenate
#   modules: losses, operations, constraints, callbacks
#


from n3fit.backends.keras_backend.internal_state import (
    set_initial_state,
    clear_backend_state,
)
from n3fit.backends.keras_backend.MetaLayer import MetaLayer
from n3fit.backends.keras_backend.MetaModel import MetaModel

from n3fit.backends.keras_backend.base_layers import (
    Input,
    Lambda,
    base_layer_selector,
    regularizer_selector,
    Concatenate,
)
from n3fit.backends.keras_backend import losses
from n3fit.backends.keras_backend import operations
from n3fit.backends.keras_backend import constraints
from n3fit.backends.keras_backend import callbacks


def select_backend(backend_name):
    """ nuke the module from orbit """
    import sys

    backends_module = sys.modules[__name__]
    # Now depenidng on the backend, we need to load the backend importer function
    if backend_name == "evolutionary_keras":
        from n3fit.backends.ga_backend.internal_state import set_backend

        set_backend(backends_module)
    elif backend_name in ["tf", "tensorflow", "keras"]:
        pass
        # from n3fit.backends.keras_backend.internal_state import set_backend
    else:
        raise ValueError(f"Backend {backend_name} not recognized")
