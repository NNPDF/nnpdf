from n3fit.backends.keras_backend import callbacks, constraints, operations
from n3fit.backends.keras_backend.MetaLayer import MetaLayer
from n3fit.backends.keras_backend.MetaModel import (
    NN_LAYER_ALL_REPLICAS,
    NN_PREFIX,
    PREPROCESSING_LAYER_ALL_REPLICAS,
    MetaModel,
)
from n3fit.backends.keras_backend.base_layers import (
    Concatenate,
    Input,
    Lambda,
    base_layer_selector,
    regularizer_selector,
)
from n3fit.backends.keras_backend.internal_state import (
    clear_backend_state,
    get_physical_gpus,
    set_eager,
    set_initial_state,
)
from n3fit.backends.keras_backend.multi_dense import MultiInitializer

print("Using Keras backend")
