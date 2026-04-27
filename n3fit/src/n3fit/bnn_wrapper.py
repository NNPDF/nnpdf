"""
Wrapper for BNN (Bayesian Neural Network) inference

This module provides utilities for:
1. Detecting if a model is a BNN (has VBDense layers)
2. Generating pseudo-replicas from BNN weight samples for plotting/analysis
"""

from n3fit.backends.keras_backend.base_layers import VBDense
from n3fit.layers.preprocessing import BayesianPreprocessing


def is_bayesian_model(pdf_model):
    """
    Check if the given pdf_model is a BNN (contains VBDense layers)

    Parameters
    ----------
    pdf_model : MetaModel
        The PDF model to check

    Returns
    -------
    bool
        True if the model is a BNN (has VBDense layers), False otherwise
    """
    vb_layers = get_vb_layers(pdf_model)
    return len(vb_layers) > 0


def _get_all_layers_recursively(container):
    """
    Recursively get all layers at any depth
    Ques.: At what layer depth of network would this hit RecursionError?
    """
    layers = [container]
    if hasattr(container, 'layers'):
        for sub_layer in container.layers:
            layers.extend(_get_all_layers_recursively(sub_layer))
    return layers


def get_vb_layers(pdf_model):
    """
    Extract all VBDense layers from a PDF model.
    Uses recursion to find VBDense at any depth in the model hierarchy.
    """
    vb_layers = []

    # Recursively get all layers at any depth
    all_layers = _get_all_layers_recursively(pdf_model)

    # Check each layer using isinstance
    for layer in all_layers:
        if isinstance(layer, VBDense):
            vb_layers.append(layer)

    return vb_layers


def get_bayesian_preprocessing(pdf_model):
    """Extract the BayesianPreprocessing layer from a PDF model, if present."""
    for layer in _get_all_layers_recursively(pdf_model):
        if isinstance(layer, BayesianPreprocessing):
            return layer
    return None


def set_model_eval(replica_model):
    """Set all VBDense and BayesianPreprocessing layers to inference mode."""
    for layer in get_vb_layers(replica_model):
        layer.eval()
    preproc = get_bayesian_preprocessing(replica_model)
    if preproc is not None:
        preproc.eval()


class BNNPredictor:
    """
    Predictor class for BNNs

    This class handles sampling from the posterior distribution
    to generate predictions with uncertainty estimates.
    """

    def __init__(self, pdf_model, n_samples=3):
        """
        Initialize the BNN predictor

        Parameters
        ----------
        pdf_model : MetaModel
            The trained PDF model with VBDense layers
        n_samples : int
            Number of samples to generate
        """
        self.pdf_model = pdf_model
        self.n_samples = n_samples
        self.vb_layers = get_vb_layers(pdf_model)
        self.bayesian_preproc = get_bayesian_preprocessing(pdf_model)

    def generate_bnn_replica(self):
        replica_models = []

        for _ in range(self.n_samples):
            replica = self.pdf_model.single_replica_generator(0)

            # Transfer VBDense posterior params
            new_vb_layers = get_vb_layers(replica)
            for parent_vb, child_vb in zip(self.vb_layers, new_vb_layers):
                child_vb.mu_w.assign(parent_vb.mu_w)
                child_vb.logsig2_w.assign(parent_vb.logsig2_w)
                child_vb.bias.assign(parent_vb.bias)

            # Transfer preprocessing alpha/beta
            replica.set_replica_weights(self.pdf_model.get_replica_weights(0), i_replica=0)

            # Fix weights and preprocessing for this replica
            set_model_eval(replica)

            replica_models.append(replica)

        return replica_models
