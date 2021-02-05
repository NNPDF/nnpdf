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


def msr_impose(fit_layer, final_pdf_layer, mode='All', scaler=None, verbose=False):
    """
    Applies the normalization due to the sum rules to the final_pdf_layer.
    The sum rules are computed in the ``fit_layer``
    while they are applied in the ``final_pdf_layer`` which is in the fk-basis.

    Parameters
    ----------
        fit_layer: layer
            The basis of the fit
        final_pdf_layer: layer
            The 14-flavours layers that goes into the fktable
        mode: str
            Which sum rules to apply (All/MSR/VSR/None)
        scaler: scaler
            Function to apply to the input. If given the input to the model
            will be a (1, None, 2) tensor where dim [:,:,0] is scaled 
    """
    # 1. Generate the fake input which will be used to integrate
    nx = int(2e3)
    xgrid, weights_array = gen_integration_input(nx)
    # 1b If a scaler is provided, scale the input xgrid
    if scaler:
        xgrid = scaler(xgrid)

    # 2. Prepare the pdf for integration
    #    for that we need to multiply several flavours with 1/x
    division_by_x = xDivide()

    def pdf_integrand(x):
        """ If a scaler is given, the division needs to take only the original input """
        x_original = operations.op_gather_keep_dims(x, -1, axis=-1)
        res = operations.op_multiply([division_by_x(x_original), fit_layer(x)])
        return res

    # 3. Now create the integration layer (the layer that will simply integrate, given some weight
    integrator = xIntegrator(weights_array, input_shape=(nx,))

    # 4. Now create the normalization by selecting the right integrations
    normalizer = MSR_Normalization(input_shape=(8,), mode=mode)

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
