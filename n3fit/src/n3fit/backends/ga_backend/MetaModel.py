"""
    MetaModel class for evolutionary_keras

    Extension of the backend MetaModel class to use evolutionary_keras.
    The MetaModel class provided by this module is a multiple inheritance
    of the keras_backend.MetaModel and the evolutionary_keras.EvolModel
"""

import n3fit.backends.keras_backend.MetaModel as tf_MetaModel
import evolutionary_keras.optimizers as Evolutionary_optimizers
from evolutionary_keras.models import EvolModel

# Add the evolutionary algorithms to the list of accepted optimizers
tf_MetaModel.optimizers["NGA"] = (
    Evolutionary_optimizers.NGA,
    {"sigma_init": 15, "population_size": 80, "mutation_rate": 0.05},
)
tf_MetaModel.optimizers["CMA"] = (
    Evolutionary_optimizers.CMA,
    {"sigma_init": 0.3, "population_size": None, "max_evaluations": None},
)

# Mix inheritance
class MetaModel(tf_MetaModel.MetaModel, EvolModel):

    accepted_optimizers = tf_MetaModel.optimizers
