"""
    The constraint module include functions to impose the momentum sum rules on the PDFs
"""
import logging
import numpy as np

from n3fit.layers import xDivide, MSR_Normalization, xIntegrator, xMultiply
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


def layer_to_integrand(layer, xgrid, eps):
    """
    Given a layer (f) and a grid in x
    return a Model corresponding to x*(f(x)-f(x-eps))

    # Arguments:
        - `layer`: layer to make into an integrand
        - `xgrid`: grid in x (numpy array)
        - `eps`: shift on x

    # Returns:
        - `modelito`: model corresponding to x*(f(x)-f(x-eps))
    """
    xgrid = xgrid.reshape(-1,1)
    # Generate the input layers
    xgrid_input = operations.numpy_to_input(xgrid, no_reshape = True)
    xgrid_input_prime = operations.numpy_to_input(xgrid - eps, no_reshape = True)
    f1 = layer(xgrid_input)
    f2 = layer(xgrid_input_prime)
    layer_diff = operations.op_subtract([f2, f1])
    # Generate a multiplier layer
    multiplier = xMultiply()
    # Generate the output layer
    output = operations.op_multiply([multiplier(xgrid_input), layer_diff])
    modelito = MetaModel([xgrid_input, xgrid_input_prime], output)
    return modelito


def compute_arclength(fitbasis_layer, nx=int(1e3)):
    """
    Given the layer with the fit basis computes the arc length

    # Arguments:
        - `fitbasis_layer`: reference to the fit basis layer
        - `nx`: number of point for the integration grid
    """
    # Generate the input layers for the xgrid and the weight
    xgrid, weights_array = gen_integration_input(nx)
    eps = xgrid[0] / 2.0
    # Generate the model to integrate
    modelito = layer_to_integrand(fitbasis_layer, xgrid, eps)
    # Compute the derivatie
    derivatives = modelito.predict(x=None, steps=1) / eps
    # Integrates the result
    derivatives_sq = pow(derivatives, 2)
    f_of_x = np.sqrt(1.0 + derivatives_sq)
    arc_lengths = np.sum(f_of_x * weights_array, axis=0)

    log.debug(
        """
        > > > Arc length:
        sigma = {0}
        g     = {1}
        v     = {2}
        v3    = {3}
        v8    = {4}""".format(
            *arc_lengths[:5]
        )
    )

    return arc_lengths
