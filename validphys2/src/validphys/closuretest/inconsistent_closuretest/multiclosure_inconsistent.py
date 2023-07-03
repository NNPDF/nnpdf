"""
closuretest/inconsistent_closuretest/multiclosure_inconsistent.py

Module containing all of the statistical estimators which are
averaged across multiple inconsistent fits. The actions
in this module are used to produce results which are plotted in
``multiclosure_inconsistent_output.py``

"""

import numpy as np

from reportengine import collect


""" To load several multiclosure fits. Useful for inconsistent closure test analysis """
multi_dataset_loader = collect("internal_multiclosure_dataset_loader", ("dataspecs",))

multi_dataset_fits_bias_replicas_variance_samples = collect(
    "dataset_fits_bias_replicas_variance_samples", ("dataspecs",)
)

multi_fits_bootstrap_dataset_bias_variance = collect(
    "fits_bootstrap_dataset_bias_variance", ("dataspecs",)
)

multi_bias_variance_resampling_dataset = collect("bias_variance_resampling_dataset", ("dataspecs",))


