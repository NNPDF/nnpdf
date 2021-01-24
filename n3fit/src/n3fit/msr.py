"""
    The constraint module include functions to impose the momentum sum rules on the PDFs
"""
import logging
import numpy as np
from scipy.interpolate import PchipInterpolator

from n3fit.layers import xDivide, MSR_Normalization, xIntegrator
from n3fit.backends import operations
from n3fit.backends import MetaModel


log = logging.getLogger(__name__)


def gen_integration_input(nx, mapping):
    """
    Generates a np.array (shaped (nx,1)) of nx elements where the
    nx/2 first elements are a logspace between 0 and 0.1
    and the rest a linspace from 0.1 to 0, which is than scaled
    using the interpolator if feature scaling is used.
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

    if mapping:
        interpolation = PchipInterpolator(mapping[0], mapping[1])
        xgrid_scaled = interpolation(np.log(xgrid.squeeze()))
        xgrid_scaled = np.expand_dims(xgrid_scaled, axis=1)

    if mapping:
        return xgrid, xgrid_scaled, weights_array
    else:
        return xgrid, weights_array


def msr_impose(fit_layer, final_pdf_layer, mapping, mode='All', verbose=False):
    """
        This function receives:
            - fit_layer: the 8-basis layer of PDF which we fit
            - final_layer: the 14-basis which is fed to the fktable
        It uses pdf_fit to compute the sum rule and returns a modified version of
        the final_pdf layer with a normalisation by which the sum rule is imposed
    """
    # 1. Generate the fake input which will be used to integrate
    nx = int(2e3)
    if mapping:
        xgrid, xgrid_scaled, weights_array = gen_integration_input(nx, mapping)
    else:
        xgrid, weights_array = gen_integration_input(nx, mapping)

    # 2. Prepare the pdf for integration
    #    for that we need to multiply several flavours with 1/x
    division_by_x = xDivide()

    def pdf_integrand(xgrid, xgrid_scaled=None):
        if mapping:
            res = operations.op_multiply([division_by_x(xgrid), fit_layer(xgrid, xgrid_scaled)])
        else:
            res = operations.op_multiply([division_by_x(xgrid), fit_layer(xgrid)])
        return res

    # 3. Now create the integration layer (the layer that will simply integrate, given some weight
    integrator = xIntegrator(weights_array, input_shape=(nx,))

    # 4. Now create the normalization by selecting the right integrations
    normalizer = MSR_Normalization(input_shape=(8,), mode=mode)

    # 5. Make the xgrid numpy array into a backend input layer so it can be given
    if mapping:
        xgrid_input_scaled = operations.numpy_to_input(xgrid_scaled)
        xgrid_input = np.expand_dims(xgrid, 0)
        xgrid_input = operations.op_convert_to_tensor(xgrid_input, dtype=xgrid_input_scaled.dtype)
        normalization = normalizer(integrator(pdf_integrand(xgrid_input, xgrid_input_scaled)))
    else:
        xgrid_input = operations.numpy_to_input(xgrid)
        normalization = normalizer(integrator(pdf_integrand(xgrid_input)))


    def ultimate_pdf(x, xgrid_scaled=None):
        if xgrid_scaled != None:
            return operations.op_multiply_dim([final_pdf_layer(x, xgrid_scaled), normalization])
        else:
            return operations.op_multiply_dim([final_pdf_layer(x), normalization])

    if verbose:
        #         only_int = integrator(pdf_integrand(xgrid_input))
        #         modelito = MetaModel(xgrid_input, only_int)
        #         result = modelito.predict(x = None, steps = 1)

        print(" > > Generating model for the inyection layer which imposes MSR")
        if mapping:
            check_integration(ultimate_pdf, xgrid_input_scaled, mapping)
        else:
            check_integration(ultimate_pdf, xgrid, mapping=None)

    # Save a reference to xgrid in ultimate_pdf, very useful for debugging
    ultimate_pdf.ref_xgrid = xgrid_input

    if mapping:
        return ultimate_pdf, xgrid_input_scaled
    else:
        return ultimate_pdf, xgrid_input


def check_integration(ultimate_pdf, integration_input, mapping):
    """
    Naive integrator for quick checks.
    Receives the final PDF layer, computes the 4 MSR and prints out the result

    Called only (for debugging purposes) by msr_impose above
    """
    nx = int(1e4)
    if mapping:
        xgrid, xgrid_scaled, weights_array = gen_integration_input(nx, mapping)
        xgrid_input_scaled = operations.numpy_to_input(xgrid_scaled)
    else:
        xgrid, weights_array = gen_integration_input(nx)
    xgrid_input = operations.numpy_to_input(xgrid)

    multiplier = xDivide(output_dim=14, div_list=range(3, 9))

    def pdf_integrand(x, x_scaled):
        res = operations.op_multiply([multiplier(x), ultimate_pdf(x, x_scaled)])
        return res

    if mapping:
        modelito = MetaModel(
            [xgrid_input, integration_input, xgrid_input_scaled],
            pdf_integrand(xgrid_input, xgrid_input_scaled),
        )
    else:
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
