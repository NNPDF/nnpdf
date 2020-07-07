"""
Penalties that can be applied to the hyperopt loss

All functions in this module have the same signature of positional arguments:

    pdf_model: :py:class:`n3fit.backends.keras_backend.MetaModel`
        model or function taking a (1, ?, 1) array as input and return a (1, ?, 14) pdf.
    stopping_object: :py:class:`n3fit.stopping.Stopping`
        object holding the information about the validation model
        and the stopping parameters

although not all penalties must use both.

And return a float to be added to the hyperscan loss.

New penalties can be added directly in this module.
The name in the runcard must match the name used in this module.
"""
import numpy as np


def saturation(pdf_model, stopping_object, n=100, min_x=1e-6, max_x=1e-4):
    """ Checks the pdf model for saturation at small x
    by checking the slope from ``min_x`` to ``max_x``.

    Parameters
    ----------
        n: int
            Number of point to evaluate the saturation
        min_x: float
            Initial point for checking the slope
        max_x: float
            Final point for checking the slope

    Example
    -------
    >>> from n3fit.hyper_optimization.penalties import saturation 
    >>> from n3fit.model_gen import pdfNN_layer_generator
    >>> fake_fl = [{'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'cbar', 's', 'sbar']]
    >>> pdf_model = pdfNN_layer_generator(nodes=[8], activations=['linear'], seed=0, flav_info=fake_fl)
    >>> isinstance(saturation(pdf_model, None), float)
    True

    """
    look_at_flavours = [1, 2]
    x = np.logspace(np.log10(min_x), np.log10(max_x), n)
    xin = x.reshape((1, -1, 1))
    y = pdf_model.predict([xin])
    extra_loss = 0.0
    for flavour in look_at_flavours:
        xpdf = y[0, :, flavour]
        slope = np.diff(xpdf) / np.diff(np.log10(x))
        pen = abs(np.mean(slope)) + np.std(slope)
        extra_loss += 1.0 / (1e-7 + pen)
    return extra_loss


def patience(pdf_model, stopping_object, alpha=1e-4):
    """ Adds a penalty for fits that have finished too soon, which
    means the number of epochs or its patience is not optimal.
    The penalty is proportional to the validation loss and will be 0
    when the best epoch is exactly at max_epoch - patience
    The ``alpha`` factor is chosen so that at 10k epochs distance
    the penalty is 2.7 * val_loss

    Parameters
    ----------
        alpha: float
            dumping factor for the exponent

    Example
    -------
    >>> from n3fit.hyper_optimization.penalties import patience
    >>> from types import SimpleNamespace
    >>> fake_stopping = SimpleNamespace(e_best_chi2=1000, stopping_patience=500, total_epochs=5000, vl_loss=2.42)
    >>> patience(None, fake_stopping, alpha=1e-4)
    3.434143467595683

    """
    epoch_best = stopping_object.e_best_chi2
    patience = stopping_object.stopping_patience
    max_epochs = stopping_object.total_epochs
    diff = abs(max_epochs - patience - epoch_best)
    vl_loss = stopping_object.vl_loss
    return vl_loss * np.exp(alpha * diff)
