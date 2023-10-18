"""
Penalties that can be applied to the hyperopt loss

Penalties in this module usually take as signature the positional arguments:

    pdf_model: :py:class:`n3fit.backends.keras_backend.MetaModel`
        model taking a ``(1, xgrid_size, 1)`` array as input
        and returning a ``(1, xgrid_size, 14, replicas)`` pdf.

    stopping_object: :py:class:`n3fit.stopping.Stopping`
        object holding the information about the validation model
        and the stopping parameters

although not all penalties use both.

And return a float to be added to the hyperscan loss.

New penalties can be added directly in this module.
The name in the runcard must match the name used in this module.
"""
import numpy as np

from n3fit.vpinterface import N3PDF, integrability_numbers
from validphys import fitveto


def saturation(pdf_model=None, n=100, min_x=1e-6, max_x=1e-4, flavors=None, **_kwargs):
    """Checks the pdf models for saturation at small x
    by checking the slope from ``min_x`` to ``max_x``.
    Sum the saturation loss of all pdf models

    Parameters
    ----------
        n: int
            Number of point to evaluate the saturation
        min_x: float
            Initial point for checking the slope
        max_x: float
            Final point for checking the slope
        flavors: list(int)
            indices of the flavors to inspect

    Example
    -------
    >>> from n3fit.hyper_optimization.penalties import saturation
    >>> from n3fit.model_gen import pdfNN_layer_generator
    >>> fake_fl = [{'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'g', 's', 'sbar']]
    >>> pdf_model = pdfNN_layer_generator(nodes=[8], activations=['linear'], seed=0, flav_info=fake_fl, fitbasis="FLAVOUR")
    >>> isinstance(saturation(pdf_model, 5), float)
    True

    """
    if flavors is None:
        flavors = [1, 2]
    x = np.logspace(np.log10(min_x), np.log10(max_x), n)
    x = np.expand_dims(x, axis=[0, -1])
    extra_loss = 0.0

    y = pdf_model.predict({"pdf_input": x})
    xpdf = y[0, :, flavors]

    delta_logx = np.diff(np.log10(x), axis=1)
    delta_xpdf = np.diff(xpdf, axis=1)
    slope = delta_xpdf / delta_logx

    pen = abs(np.mean(slope, axis=1)) + np.std(slope, axis=1)

    # sum over flavors
    # Add a small offset to avoid ZeroDivisionError
    extra_loss += np.sum(1.0 / (1e-7 + pen), axis=0)
    return extra_loss


def patience(stopping_object=None, alpha=1e-4, **_kwargs):
    """Adds a penalty for fits that have finished too soon, which
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
    >>> patience(fake_stopping, alpha=1e-4)
    3.434143467595683

    """
    epoch_best = np.take(stopping_object.e_best_chi2, 0)
    patience = stopping_object.stopping_patience
    max_epochs = stopping_object.total_epochs
    diff = abs(max_epochs - patience - epoch_best)
    vl_loss = np.take(stopping_object.vl_chi2, 0)
    return vl_loss * np.exp(alpha * diff)


def integrability(pdf_model=None, **_kwargs):
    """Adds a penalty proportional to the value of the integrability integration
    It adds a 0-penalty when the value of the integrability is equal or less than the value
    of the threshold defined in validphys::fitveto

    The penalty increases exponentially with the growth of the integrability number

    Example
    -------
    >>> from n3fit.hyper_optimization.penalties import integrability
    >>> from n3fit.model_gen import pdfNN_layer_generator
    >>> fake_fl = [{'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'g', 's', 'sbar']]
    >>> pdf_model = pdfNN_layer_generator(nodes=[8], activations=['linear'], seed=0, flav_info=fake_fl, fitbasis="FLAVOUR")
    >>> isinstance(integrability(pdf_model), float)
    True

    """
    pdf_instance = N3PDF(pdf_model.split_replicas())
    integ_values = integrability_numbers(pdf_instance)
    integ_overflow = np.sum(integ_values[integ_values > fitveto.INTEG_THRESHOLD])
    if integ_overflow > 50.0:
        # before reaching an overflow, just give a stupidly big number
        return np.exp(50.0)
    return np.exp(integ_overflow) - 1.0
