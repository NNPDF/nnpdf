"""
    The constraint module include functions to impose the momentum sum rules on the PDFs
"""
import logging
import numpy as np

from validphys.arclength import arc_length_core_computation

from n3fit.layers import xDivide, MSR_Normalization, xIntegrator
from n3fit.backends import operations
from n3fit.backends import MetaModel


log = logging.getLogger(__name__)


def gen_integration_input(nx):
    """
    Generates a np.array (shaped (nx,1)) of nx elements where the
    nx/2 first elements are a logspace between 0 and 0.1
    and the rest a linspace from 0.1 to 0
    """
    lognx = int(nx / 2)
    linnx = int(nx - lognx)
    xgrid_log = np.logspace(-9, -1, lognx + 1)
    xgrid_lin = np.linspace(0.1, 1, linnx)
    xgrid = np.concatenate([xgrid_log[:-1], xgrid_lin]).reshape(nx, 1)

    spacing = [0.0]
    for i in range(1, nx):
        spacing.append(np.abs(xgrid[i - 1] - xgrid[i]))
    spacing.append(0.0)

    weights = []
    for i in range(nx):
        weights.append((spacing[i] + spacing[i + 1]) / 2.0)
    weights_array = np.array(weights).reshape(nx, 1)

    return xgrid, weights_array


def msr_impose(fit_layer, final_pdf_layer):
    """
        This function receives:
            - fit_layer: the 8-basis layer of PDF which we fit
            - final_layer: the 14-basis which is fed to the fktable
        It uses pdf_fit to compute the sum rule and returns a modified version of
        the final_pdf layer with a normalisation by which the sum rule is imposed
    """
    # 1. Generate the fake input which will be used to integrate
    nx = int(2e3)
    xgrid, weights_array = gen_integration_input(nx)

    # 2. Prepare the pdf for integration
    #    for that we need to multiply several flavours with 1/x
    division_by_x = xDivide()

    def pdf_integrand(x):
        res = operations.op_multiply([division_by_x(x), fit_layer(x)])
        return res

    # 3. Now create the integration layer (the layer that will simply integrate, given some weight
    integrator = xIntegrator(weights_array, input_shape=(nx,))

    # 4. Now create the normalization by selecting the right integrations
    normalizer = MSR_Normalization(input_shape=(8,))

    # 5. Make the xgrid numpy array into a backend input layer so it can be given
    xgrid_input = operations.numpy_to_input(xgrid)
    normalization = normalizer(integrator(pdf_integrand(xgrid_input)))

    def ultimate_pdf(x):
        return operations.op_multiply_dim([final_pdf_layer(x), normalization])

    # Save a reference to xgrid in ultimate_pdf, very useful for debugging
    ultimate_pdf.ref_xgrid = xgrid_input

    return ultimate_pdf, xgrid_input


def compute_arclength(pdf_function):
    """
    Receives a PDF function and returns an array with
    the computed arclenghts.
    The input function can receive an array of x and return
    an array of f number of flavours per x
    Parameters
    ----------
        pdf_function: function
            A function that given a value on (x) returns a value of x*pdf(x)
            per flavour
    """

    def pdf_values(qarr, flarr, xgrid):
        """ Generate a valid validphys grid_value function
        n3fit fits produce PDFs at a fixed Q in the evolution basis.
        The Q is ignored.
        The xgrid variable instead is a tuple where the second element
        is the x in which to compute the PDF values.
        """
        ret = pdf_function(np.expand_dims(xgrid[1], -1))
        # Select desired flavours
        ret = ret.T[flarr]
        # Ensure the output is [replicas][flavours][x][Q]
        return np.expand_dims(ret, (0, -1))

    qignore = 42
    # TODO
    # n3fit returns the PDF in the flavour order required by the fktables
    # so there should be somewhere a dictionary defining that so that we can just
    # say here the flavours we want by name
    flavour_selection = [1, 2, 3, 4, 5, 9, 10, 11]
    res = arc_length_core_computation(pdf_values, qignore, flavour_selection)
    return res[0]
