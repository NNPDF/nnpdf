"""
Penalties that can be applied to the hyperopt loss

All functions in this module will receive a `pdf` model which takes a (1, ?, 1) array as input
and return a (1, ?, 14) pdf.
All functions return a float to be added to the loss
"""
import numpy as np

def saturation(pdf_model, *args, n = 100, min_x = 1e-6, max_x = 1e-4):
    """ Checks the pdf model for saturation at small x
    by checking the slope from `min_x` to `max_x`
    """
    look_at_flavours = [1, 2]
    x = np.logspace(np.log10(min_x), np.log10(max_x), n)
    xin = x.reshape( (1,-1,1) )
    y = pdf_model.predict([xin])
    extra_loss = 0.0
    for flavour in look_at_flavours:
        xpdf = y[0,:,flavour]
        slope = np.diff(xpdf)/np.diff(np.log10(x))
        pen = abs(np.mean(slope)) + np.std(slope)
        extra_loss += 1.0/(1e-7 + pen)
    return extra_loss

def patience(pdf_model, stopping_object, alpha = 1e-4):
    """ Adds a penalty for fits that have finished too soon, which
    means the number of epochs or its patience is not optimal.
    The penalty is proportional to the validation loss and will be 0
    when the best epoch is exactly at max_epoch - patience
    The alpha factor is chosen so that at 10k epochs distance
    the penalty is 2.7 * val_loss
    """
    epoch_best = stopping_object.e_best_chi2,
    patience = stopping_object.stopping_patience
    max_epochs = stopping_object.total_epochs
    diff = abs(max_epochs - patience - epoch_best)
    vl_loss = stopping_object.vl_loss
    return vl_loss*np.exp(alpha*diff)

