"""
    The constraint module include functions to impose the momentum sum rules on the PDFs
"""
import logging
import numpy as np

from n3fit.layers import xDivide, MSR_Normalization, xIntegrator
from n3fit.backends import operations
from n3fit.backends import MetaModel
from validphys.arclength import arc_lengths, integrability_number


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


def msr_impose(fit_layer, final_pdf_layer, verbose=False):
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

    if verbose:
        #         only_int = integrator(pdf_integrand(xgrid_input))
        #         modelito = MetaModel(xgrid_input, only_int)
        #         result = modelito.predict(x = None, steps = 1)

        print(" > > Generating model for the inyection layer which imposes MSR")
        check_integration(ultimate_pdf, xgrid_input)

    # Save a reference to xgrid in ultimate_pdf, very useful for debugging
    ultimate_pdf.ref_xgrid = xgrid_input

    return ultimate_pdf, xgrid_input


def check_integration(ultimate_pdf, integration_input):
    """
    Naive integrator for quick checks.
    Receives the final PDF layer, computes the 4 MSR and prints out the result

    Called only (for debugging purposes) by msr_impose above
    """
    nx = int(1e4)
    xgrid, weights_array = gen_integration_input(nx)
    xgrid_input = operations.numpy_to_input(xgrid)

    multiplier = xDivide(output_dim=14, div_list=range(3, 9))

    def pdf_integrand(x):
        res = operations.op_multiply([multiplier(x), ultimate_pdf(x)])
        return res

    modelito = MetaModel([xgrid_input, integration_input], pdf_integrand(xgrid_input))
    modelito.summary()
    result = modelito.predict(x=None, steps=1)

    result_weighted = result * weights_array
    result_integrated = np.sum(result_weighted, axis=0)

    msr = result_integrated[1] + result_integrated[2]
    v = result_integrated[3]
    v3 = result_integrated[4]
    v8 = result_integrated[5]
    print(
        """
     > > > Int from 0 to 1 of:
    x*g(x) + x*sigma(x) = {0}
    v                   = {1}
    v3                  = {2}
    v8                  = {3}""".format(
            msr, v, v3, v8
        )
    )


# TODO: once the writer is absorbed into vp this function can dissapear as well
def compute_arclength(n3pdf):
    """
    Given the layer with the fit basis computes the arc length

    Parameters
    ----------
        pdf_function: function
            pdf function has received by the writer or ``pdf_model``

    Example
    -------

    >>> from n3fit.vpinterface import N3PDF
    >>> from n3fit.model_gen import pdfNN_layer_generator
    >>> from n3fit.msr import compute_arclength
    >>> fake_fl = [{'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'cbar', 's', 'sbar']]
    >>> pdf_model = pdfNN_layer_generator(nodes=[8], activations=['linear'], seed=0, flav_info=fake_fl)
    >>> n3pdf = N3PDF(pdf_model)
    >>> res = compute_arclength(n3pdf)
    """
    ret = arc_lengths(n3pdf, [1.65], "evolution", ["sigma", "gluon", "V", "V3", "V8"])
    return ret.stats.central_value()

def compute_integrability_number(n3pdf):
    """
    Given the layer with the fit basis computes the integrability number
    for the distributions V, V3, V8, T3, T8

    Parameters
    ----------
        n3pdf: pdf function has received by the writer or ``pdf_model``
    """
    ret = integrability_number(n3pdf, [1.65], flavours=['V', 'T3', 'V3', 'T8', 'V8'])
    return ret