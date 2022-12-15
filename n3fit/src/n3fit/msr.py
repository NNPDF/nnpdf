"""
    The constraint module include functions to impose the momentum sum rules on the PDFs
"""
import logging
import numpy as np

from n3fit.layers import xDivide, MSR_Normalization, xIntegrator
from n3fit.backends import operations as op


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


def msr_impose(nx=int(2e3), mode='All', scaler=None, photon=None):
    """
        This function receives:
        Generates a function that applies a normalization layer to the fit.
            - fit_layer: the 8-basis layer of PDF which we fit
        The normalization is computed from the direct output of the NN (so the 7,8-flavours basis)
            - final_layer: the 14-basis which is fed to the fktable
        and it is applied to the input of the fktable (i.e., to the 14-flavours fk-basis).
        It uses pdf_fit to compute the sum rule and returns a modified version of
        the final_pdf layer with a normalisation by which the sum rule is imposed

        Parameters
        ----------
            nx: int
                number of points for the integration grid, default: 2000
            mode: str
                what sum rules to compute (MSR, VSR or All), default: All
            scaler: scaler
                Function to apply to the input. If given the input to the model
                will be a (1, None, 2) tensor where dim [:,:,0] is scaled 
    """

    # 1. Generate the fake input which will be used to integrate
    xgrid, weights_array = gen_integration_input(nx)
    # 1b If a scaler is provided, scale the input xgrid
    if scaler:
        xgrid = scaler(xgrid)

    # 2. Prepare the pdf for integration
    #    for that we need to multiply several flavours with 1/x
    division_by_x = xDivide()
    # 3. Now create the integration layer (the layer that will simply integrate, given some weight
    integrator = xIntegrator(weights_array, input_shape=(nx,))

    # 3.1 If a photon is given, compute the photon component of the MSR
    photon_c = 0.0
    if photon is not None:
        photon_c = photon.integrate()

    # 4. Now create the normalization by selecting the right integrations
    normalizer = MSR_Normalization(mode=mode, photon_contribution=photon_c)

    # 5. Make the xgrid array into a backend input layer so it can be given to the normalization
    xgrid_input = op.numpy_to_input(xgrid, name="integration_grid")
    # Finally prepare a function which will take as input the output of the PDF model
    # and will return it appropiately normalized.
    def apply_normalization(layer_pdf):
        """
            layer_pdf: output of the PDF, unnormalized, ready for the fktable
        """
        x_original = op.op_gather_keep_dims(xgrid_input, -1, axis=-1)
        pdf_integrand = op.op_multiply([division_by_x(x_original), layer_pdf(xgrid_input)])
        normalization = normalizer(integrator(pdf_integrand))

        def ultimate_pdf(x):
            return layer_pdf(x)*normalization

        return ultimate_pdf

    return apply_normalization, xgrid_input
