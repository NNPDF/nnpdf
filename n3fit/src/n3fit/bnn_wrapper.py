"""
Wrapper for BNN (Bayesian Neural Network) inference 

This module provides utilities for:
1. Detecting if a model is a BNN (has VBDense layers)
2. Generating pseudo-replicas from BNN weight samples for plotting/analysis
"""

import numpy as np
import tensorflow as tf


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


def get_vb_layers(pdf_model):
    """
    Extract all VBDense layers from a PDF model.
    
    Parameters
    ----------
    pdf_model : MetaModel
        The PDF model to extract layers from
        
    Returns
    -------
    list
        List of VBDense layer instances
    """
    vb_layers = []
    
    # Check all layers in the model
    for layer in pdf_model.layers:
        if hasattr(layer, '__class__') and layer.__class__.__name__ == 'VBDense':
            vb_layers.append(layer)
        
        # Also check nested layers (for MetaModel structure)
        if hasattr(layer, 'layers'):
            for sub_layer in layer.layers:
                if hasattr(sub_layer, '__class__') and sub_layer.__class__.__name__ == 'VBDense':
                    vb_layers.append(sub_layer)
    
    return vb_layers


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
        
    def reset_random(self):
        """Reset the random state for each VBDense layer."""
        for vb_layer in self.vb_layers:
            vb_layer.reset_random()
            
    def generate_predictions(self, xinput):
        """
        Generate predictions via sampling
        
        This uses the same machinery as model_trainer.py (_model_generation method)
        to integrate with n3fit's observable layers.
        
        Parameters
        ----------
        xinput : InputInfo
            The input info containing the pdf_input layer and split layer
            
        Returns
        -------
        predictions : np.ndarray
            Array of shape (n_samples, batch_size, n_replicas, n_xgrid, n_flavours)
            containing all replica predictions
        """
        predictions = []
        
        for i in range(self.n_samples):
            # Reset random weights for each VBDense layer
            self.reset_random()
            
            # Get PDF using apply_as_layer (same pattern as model_trainer)
            full_model_input_dict, full_pdf = self.pdf_model.apply_as_layer(
                {"pdf_input": xinput.input}
            )
            
            # Evaluate the PDF
            pdf_output = full_pdf.numpy()
            predictions.append(pdf_output)
            
        predictions = np.array(predictions)
        return predictions
    
    def generate_predictions_simple(self, x_test):
        """
        Simplified prediction using model directly.
        
        Parameters
        ----------
        x_test : tf.Tensor or np.ndarray
            Input x values with shape (1, n_x, 1)
            
        Returns
        -------
        predictions : np.ndarray
            Array of replica predictions
        """
        predictions = []
        
        for i in range(self.n_samples):
            # Reset random weights
            self.reset_random()
            
            # Forward pass with training=False for inference mode
            pred = self.pdf_model(x_test, training=False)
            predictions.append(pred.numpy())
            
        predictions = np.array(predictions)
        return predictions
    
def bayesian_inference_with_observables(model_trainer, pdf_model, n_samples=3):
    """
    Full pipeline: PDF inference -> observable computation -> chi2
    
    This function demonstrates how to use the BNN with n3fit's observable
    framework, similar to how model_trainer.py works
    
    Parameters
    ----------
    model_trainer : ModelTrainer
        The trained ModelTrainer containing the experimental setup
    pdf_model : MetaModel
        The trained PDF model with VBDense layers
    n_samples : int
        Number of replica samples to be drawn
        
    Returns
    -------
    mean_pdf: np.ndarray
            Mean PDF predictions across all samples
    epistemic: np.ndarray
            Epistemic uncertainty (std of predictions)
    """
    from n3fit.bnn_wrapper import BNNPredictor
    
    # Generate input grid
    xinput = model_trainer._xgrid_generation()
    
    # Create predictor
    bnn_predictor = BNNPredictor(pdf_model, n_samples=n_samples)
    
    # Get predictions
    predictions = bnn_predictor.generate_predictions_simple(xinput)
    
    # Compute statistics
    mean_pred = np.mean(predictions, axis=0)
    epistemic = np.std(predictions, axis=0)
    
    return mean_pred, epistemic

