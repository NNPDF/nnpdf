"""
    The constraint module include functions to impose the momentum sum rules on the PDFs
"""
import logging
import numpy as np

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


def msr_impose(nx=int(2e3), basis_size=8, mode='All'):
    """
        This function receives:
            - fit_layer: the 8-basis layer of PDF which we fit
            - final_layer: the 14-basis which is fed to the fktable
        It uses pdf_fit to compute the sum rule and returns a modified version of
        the final_pdf layer with a normalisation by which the sum rule is imposed
    """
    # 1. Generate the fake input which will be used to integrate
    xgrid, weights_array = gen_integration_input(nx)

    # 2. Prepare the pdf for integration
    #    for that we need to multiply several flavours with 1/x
    division_by_x = xDivide()

    # 3. Now create the integration layer (the layer that will simply integrate, given some weight
    integrator = xIntegrator(weights_array, input_shape=(nx,))

    # 4. Now create the normalization by selecting the right integrations
    normalizer = MSR_Normalization(input_shape=(basis_size,), mode=mode)

    # 5. Make the xgrid array into a backend input layer so it can be given to the normalization
    xgrid_input = operations.numpy_to_input(xgrid)

    # Now parepare a function that takes as input the 8-flavours output of the NN
    # and the 14-flavours after the fk rotation and returns a 14-flavours normalized output
    # note + TODO:
    # the idea was that the normalization should always be applied at the fktable 14-flavours
    # and always computed at the output of the NN (in case one would like to compute it differently)
    # don't think it is a good idea anymore and should be changed to act only on the output to the fktable
    # but will be dealt with in the future.
                            # fitlayer        #final_pdf
    def apply_normalization(layer_fitbasis, layer_pdf):
        """
            layer_fitbasis: output of the NN
            layer_pdf: output for the fktable
        """

        pdf_integrand = operations.op_multiply([division_by_x(xgrid_input), layer_fitbasis(xgrid_input)])
        normalization = normalizer(integrator(pdf_integrand))

        def ultimate_pdf(x):
            return operations.op_multiply([layer_pdf(x), normalization])

        return ultimate_pdf

    return apply_normalization, xgrid_input


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
