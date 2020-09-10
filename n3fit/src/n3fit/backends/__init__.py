from n3fit.backends.keras_backend.internal_state import (
    set_initial_state,
    clear_backend_state,
)
from n3fit.backends.keras_backend.MetaLayer import MetaLayer

from n3fit.backends.keras_backend.base_layers import (
    Input,
    concatenate,
    Lambda,
    base_layer_selector,
    regularizer_selector,
    Concatenate,
)
from n3fit.backends.keras_backend import losses
from n3fit.backends.keras_backend import operations
from n3fit.backends.keras_backend import constraints
from n3fit.backends.keras_backend import callbacks

# Don't import the Model until it needs to be imported
class _MetaModel:
    def __init__(self, backend="tensorflow"):
        self.backend = backend
        from n3fit.backends.keras_backend.MetaModel import MetaModel

        self.meta_model = MetaModel

    def enable_ga(self):
        try:
            from n3fit.backends.ga_backend.MetaModel import MetaModel

            self.meta_model = MetaModel
            self.backend = "evolutionary_keras"
        except ModuleNotFoundError:
            raise ModuleNotFoundError("Install `evolutionary_keras` to use this backend")

    def __call__(self, *args, **kwargs):
        return self.meta_model(*args, **kwargs)


MetaModel = _MetaModel()
