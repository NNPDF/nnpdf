from validphys.api import API

import numpy as np

# mpmath for algebraic arbitrary precision
from mpmath import gamma
from scipy.optimize import curve_fit


def set_template(fit, t0pdfset):
    # didn't find a more practical way
    template = {
        "theory": {"from_": "fit"},
        "theoryid": {"from_": "theory"},
        "use_cuts": "fromfit",
        "experiments": {"from_": "fit"},
        "fit": fit,
        "use_t0": True,
        "t0pdfset": t0pdfset,
    }

    return template


# 'pdf' is the converted Hessian set of 'fit' used in 'template'
def delta_chi2_eigs(pdf, template):
    """
    Return cv and std dev of delta_chi2 over all eigenvectors.
    """
    # delta_chi2_hessian implemented in validphys/dataplots.py
    return API.delta_chi2_hessian(pdf=pdf, **template)


def _G(m, N):
    """
    Compute the function G_n(m) which appears in the coefficients. Store the computed
    values to calculate them only one single time.
    """
    name = str(m) + str(N)
    if name in G_values:
        return G_values[name]
    else:
        # gamma function with arbitrary precision from mpmath
        result = gamma((N - 1) / 2 + m / 2) / gamma((N - 1) / 2) * np.power(2 / (N - 1), m / 2)
        # gamma return an mpf type
        G_values[name] = np.array(result).astype(np.float64)
        return G_values[name]


G_values = {}  # empty dictionary to store values of the function G_n(m)


def coeff_n_indep(a_n, b_n, c_n, d_n, nrep, sign=None):
    """
    Return the coefficients a, b, c, d nrep independent.

    Coefficients a_n and c_n has an ambiguity of sign.
    """
    a = a_n / _G(1.0, nrep)
    if sign is not None:
        a *= -1

    b = b_n

    c = c_n / (_G(3.0, nrep) + 3 / nrep * _G(1.0, nrep))
    if sign is not None:
        c *= -1

    d = d_n * (nrep * (nrep - 1)) / (nrep ** 2 + 7 * nrep - 6)

    return [a, b, c, d]


def coeff_n_dep(a, b, c, d, nrep, sign=None):
    """
    Return the coefficients a_n, b_n, c_n, d_n nrep dependent.
    """
    a_n = a * _G(1.0, nrep)
    if sign is not None:
        a_n *= -1

    b_n = b

    c_n = c * (_G(3.0, nrep) + 3 / nrep * _G(1.0, nrep))
    if sign is not None:
        c_n *= -1

    d_n = d * (nrep ** 2 + 7 * nrep - 6) / (nrep * (nrep - 1))

    return [a_n, b_n, c_n, d_n]


def model_std(a, b, c, d, nrep):
    """
    Return the std dev of the coefficients nrep dependent as described by the model.
    """
    std_a_n = a * np.sqrt(1 - _G(1.0, nrep) ** 2)

    std_b_n = b * np.sqrt(2 * (3 * nrep - 2) / (nrep * (nrep - 1)))

    std_c_n = c * np.sqrt(
        (nrep ** 3 + 19 * nrep ** 2 + 3 * nrep - 15) / (nrep * (nrep - 1) ** 2)
        - (_G(3.0, nrep) + 3 / nrep * _G(1.0, nrep)) ** 2
    )

    std_d_n = d * np.sqrt(
        (nrep + 1) * (nrep + 3) * (nrep + 5) / (nrep - 1) ** 3
        + 28 * (nrep + 1) * (nrep + 3) / (nrep * (nrep - 1) ** 2)
        + 140 * (nrep + 1) / (nrep ** 2 * (nrep - 1))
        + 540 / nrep ** 3
        - ((nrep ** 2 + 7 * nrep - 6) / (nrep * (nrep - 1))) ** 2
    )

    return [std_a_n, std_b_n, std_c_n, std_d_n]


def model_delta_chi2(x, a_n, b_n, c_n, d_n):

    return a_n * x + b_n * x ** 2 + c_n * x ** 3 + d_n * x ** 4


def fit_delta_chi2(x, y, s):
    """
    Implementation of delta_chi2 fit.

    Return the fitted coefficients nrep dependent, and their std dev estimated by the fit.
    """
    # popt are the estimated coeffs a_n, b_n, c_n, d_n
    popt, pcov = curve_fit(model_delta_chi2, x, y, sigma=s)

    return popt, np.sqrt(np.diag(pcov))
