
import numpy as np

def symmetrize_errors(delta_plus, delta_minus):
    r"""Compute the symmetrized uncertainty and the shift in data point.

    Parameters
    ----------
    delta_plus : float
        The top/plus uncertainty with sign
    delta_minus : float
        The bottom/minus uncertainty with sign

    Returns
    -------
    se_delta : float
        The value to be added to the data point
    se_sigma : float
        The symmetrized uncertainty to be used in commondata

    """
    semi_diff = (delta_plus + delta_minus) / 2
    average = (delta_plus - delta_minus) / 2
    se_delta = semi_diff
    se_sigma = np.sqrt(average * average + 2 * semi_diff * semi_diff)
    return se_delta, se_sigma

