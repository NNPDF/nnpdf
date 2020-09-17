# Make the chosen backend available from n3fit.backends
# These includes:
#   set_initial_state, clear_backend_state
#   meta classes: MetaLayer, MetaModel
#   modules: backend_layers, losses, operations, constraints, callbacks
#

# If no backend has been set the default is the tensorflow backend
from sys import modules as _modules
_backends_module = _modules[__name__]

from n3fit.backends.keras_backend.internal_state import set_backend as set_tf
set_tf(_backends_module)

def select_backend(backend_name):
    """ Select the appropiate backend by overriding this module 
    Default is understood as TensorFlow
    """
    if backend_name is None:
        backend_name = "tf"
    try:
        backend_name = backend_name.lower()
    except AttributeError:
        raise ValueError(f"select_backend accepts only strings, received: {backend_name}")
    # Now depenidng on the backend, we need to load the backend importer function
    if backend_name == "evolutionary_keras":
        from n3fit.backends.ga_backend.internal_state import set_backend

        set_backend(_backends_module)
    elif backend_name in ["tf", "tensorflow", "keras"]:
        set_tf(_backends_module)
    else:
        raise ValueError(f"Backend {backend_name} not recognized")
