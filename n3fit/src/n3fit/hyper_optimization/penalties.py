"""
Penalties that can be applied to the hyperopt loss

Penalties in this module usually take as signature the positional arguments:

    pdf_models: list(:py:class:`n3fit.backends.keras_backend.MetaModel`)
        list of models or functions taking a ``(1, xgrid_size, 1)`` array as input
        and returns a ``(1, xgrid_size, 14)`` pdf.

    stopping_object: :py:class:`n3fit.stopping.Stopping`
        object holding the information about the validation model
        and the stopping parameters

although not all penalties use both.

And return a float to be added to the hyperscan loss.

New penalties can be added directly in this module.
The name in the runcard must match the name used in this module.
"""
from typing import List, Optional

import numpy as np
from numpy.typing import NDArray

from n3fit.vpinterface import N3PDF, integrability_numbers
from validphys import fitveto


def saturation(
    pdf_models: List,
    n: int = 100,
    min_x: float = 1e-6,
    max_x: float = 1e-4,
    flavors: Optional[List[int]] = None,
    **_kwargs
) -> NDArray:
    """
    Checks the pdf models for saturation at small x
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

    Returns
    -------
        NDArray
            array of saturation penalties for each replica

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
    xin = np.expand_dims(x, axis=[0, -1])
    extra_losses = []
    for pdf_model in pdf_models:
        y = pdf_model.predict({"pdf_input": xin})
        xpdf = y[0, :, flavors]
        slope = np.diff(xpdf) / np.diff(np.log10(x))
        pen = abs(np.mean(slope, axis=1)) + np.std(slope, axis=1)
        # Add a small offset to avoid ZeroDivisionError
        extra_losses.append(np.sum(1.0 / (1e-7 + pen)))
    return np.array(extra_losses)


def patience(stopping_object, alpha: float = 1e-4, **_kwargs) -> NDArray:
    """
    Adds a penalty for fits that have finished too soon, which
    means the number of epochs or its patience is not optimal.
    The penalty is proportional to the validation loss and will be 0
    when the best epoch is exactly at max_epoch - patience
    The ``alpha`` factor is chosen so that at 10k epochs distance
    the penalty is 2.7 * val_loss

    Parameters
    ----------
        alpha: float
            dumping factor for the exponent

    Returns
    -------
        NDArray
            patience penalty for each replica

    Example
    -------
    >>> from n3fit.hyper_optimization.penalties import patience
    >>> from types import SimpleNamespace
    >>> fake_stopping = SimpleNamespace(e_best_chi2=1000, stopping_patience=500, total_epochs=5000, vl_loss=2.42)
    >>> patience(fake_stopping, alpha=1e-4)
    3.434143467595683

    """
    epoch_best = np.array(stopping_object.e_best_chi2)
    patience = stopping_object.stopping_patience
    max_epochs = stopping_object.total_epochs
    diff = abs(max_epochs - patience - epoch_best)
    vl_loss = np.array(stopping_object.vl_chi2)
    return vl_loss * np.exp(alpha * diff)


def integrability(pdf_models: List, **_kwargs) -> NDArray:
    """
    Adds a penalty proportional to the value of the integrability integration
    It adds a 0-penalty when the value of the integrability is equal or less than the value
    of the threshold defined in validphys::fitveto

    The penalty increases exponentially with the growth of the integrability number

    Returns
    -------
        NDArray
            array of integrability penalties for each replica

    Example
    -------
    >>> from n3fit.hyper_optimization.penalties import integrability
    >>> from n3fit.model_gen import pdfNN_layer_generator
    >>> fake_fl = [{'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'g', 's', 'sbar']]
    >>> pdf_model = pdfNN_layer_generator(nodes=[8], activations=['linear'], seed=0, flav_info=fake_fl, fitbasis="FLAVOUR")
    >>> isinstance(integrability(pdf_model), float)
    True

    """
    pdf_instance = N3PDF(pdf_models)
    integ_values = integrability_numbers(pdf_instance)

    # set components under the threshold to 0
    integ_values[integ_values <= fitveto.INTEG_THRESHOLD] = 0.0

    # sum over flavours
    integ_overflow = np.sum(integ_values, axis=1)

    # limit components to 50 to avoid overflow
    integ_overflow[integ_overflow > 50.0] = 50.0

    return np.exp(integ_overflow) - 1.0
