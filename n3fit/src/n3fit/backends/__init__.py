from n3fit.backends.keras_backend.internal_state import (
    set_initial_state,
    clear_backend_state,
    set_eager
)
from n3fit.backends.keras_backend.MetaLayer import MetaLayer
from n3fit.backends.keras_backend.MetaModel import MetaModel
from n3fit.backends.keras_backend.base_layers import (
    Input,
    concatenate,
    Lambda,
    base_layer_selector,
    regularizer_selector,
    Concatenate,
)
from n3fit.backends.keras_backend import operations
from n3fit.backends.keras_backend import constraints
from n3fit.backends.keras_backend import callbacks

print("Using Keras backend")
