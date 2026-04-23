"""
Wrapper for BNN (Bayesian Neural Network) inference

This module provides utilities for:
1. Detecting if a model is a BNN (has VBDense layers)
2. Generating pseudo-replicas from BNN weight samples for plotting/analysis
"""

import numpy as np
import tensorflow as tf
from n3fit.backends.keras_backend.base_layers import VBDense

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

class BNNPredictor:
    """
    Predictor class for BNNs

    This class handles sampling from the posterior distribution
    to generate predictions with uncertainty estimates.
    """

    def __init__(self, pdf_model, training_model, n_samples=3):
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
        self.training_model = training_model

    def reset_random(self):
        """Reset the random state for each VBDense layer."""
        for vb_layer in self.vb_layers:
            vb_layer.reset_random()

    def eval(self):
        """Evaluate the model in inference mode (training=False)."""
        for vb_layer in self.vb_layers:
            vb_layer.eval()

    def train(self):
        """Set the model to training mode (training=True)."""
        for vb_layer in self.vb_layers:
            vb_layer.train()

    def generate_bnn_replica_old(self):
        replica_models =[]
        for i in range(self.n_samples):
            self.reset_random()
            self.eval()

            replica = self.pdf_model.single_replica_generator(0)
            #replica.set_replica_weights(self.pdf_model.get_replica_weights(0), i_replica=0)
            replica_models.append(replica)

        return replica_models
    
    def generate_bnn_replica(self):
        replica_models = []
        x_input = {"pdf_input": self.training_model.x_in["pdf_input"]}
        
        for i in range(self.n_samples):
            # Create a fresh single-replica model with the same architecture
            replica = self.pdf_model.single_replica_generator(0)
            
            # Transfer trained posterior parameters (mu, rho) from parent VBDense layers
            new_vb_layers = get_vb_layers(replica)
            for parent_vb, child_vb in zip(self.vb_layers, new_vb_layers):
                child_vb.mu_w.assign(parent_vb.mu_w)
                child_vb.logsig2_w.assign(parent_vb.logsig2_w)
                child_vb.bias.assign(parent_vb.bias)
            
            # Also transfer preprocessing weights so the replica is fully consistent
            replica.set_replica_weights(
                self.pdf_model.get_replica_weights(0), i_replica=0
            )
            
            replica_models.append(replica)
        
        for i, replica in enumerate(replica_models[:3]):
            vb = get_vb_layers(replica)[0]
            print(f"Replica {i}: mu_w mean = {vb.mu_w.numpy().mean():.6f}")
        
        return replica_models

    def generate_bnn_replica_return_nparray(self):
        # returns numpy arrays (output of predict()), 
        # but N3PDF expects MetaModel objects same as
        # what split_replicas() returns in the non-BNN path.
        # BUT THIS IS THE BASE LOGIC FOR REPLICA SAMPLING
        
        replica_outputs = []
        # Extract x-grid form the training model's compiled inputs
        x_input = {"pdf_input": self.training_model.x_in["pdf_input"]}
        for i in range(self.n_samples):
            # Each call draws fresh epsilon ~ N(0,I), giving a new weight sample
            # from the trained posterior N(mu_w, exp(0.5*logsig2_w))
            output = self.pdf_model.predict(x=x_input)
            replica_outputs.append(output)
        return replica_outputs