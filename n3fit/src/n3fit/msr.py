"""
    The constraint module include functions to impose the momentum sum rules on the PDFs
"""
import logging
import numpy as np

from n3fit.layers import xDivide, MSR_Normalization, xIntegrator
from n3fit.backends import operations as op

from n3fit.backends import MetaModel, Lambda, Input


log = logging.getLogger(__name__)


def generate_msr_model_and_grid(
        output_dim: int = 14,
        mode: str = "ALL",
        nx: int = int(2e3),
        scaler=None,
        **kwargs) -> MetaModel:
    """
    Generates a model that applies the sum rules to the PDF.

    Parameters
    ----------
    output_dim: int
        Number of flavours of the output PDF
    mode: str
        Mode of sum rules to apply. It can be:
            - "ALL": applies both the momentum and valence sum rules
            - "MSR": applies only the momentum sum rule
            - "VSR": applies only the valence sum rule
    nx: int
        Number of points of the integration grid
    scaler: Scaler
        Scaler to be applied to the PDF before applying the sum rules

    Returns
    -------
    model: MetaModel
        Model that applies the sum rules to the PDF
        It takes as inputs:
            - pdf_x: the PDF output of the model
            - pdf_xgrid_integration: the PDF output of the model evaluated at the integration grid
            - xgrid_integration: the integration grid
        It returns the PDF with the sum rules applied
    xgrid_integration: dict
        Dictionary with the integration grid, with:
            - values: the integration grid
            - input: the input layer of the integration grid
    """
    # 0. Prepare input layers to MSR model
    pdf_x = Input(shape=(None, output_dim), batch_size=1, name="pdf_x")
    pdf_xgrid_integration = Input(shape=(nx, output_dim), batch_size=1, name="pdf_xgrid_integration")

    # 1. Generate the grid and weights that will be used to integrate
    xgrid_integration, weights_array = gen_integration_input(nx)
    # 1b If a scaler is provided, scale the input xgrid
    if scaler:
        xgrid_integration = scaler(xgrid_integration)

    # Turn into input layer. Specify a custom shape here using None,
    # so shapes will display properly in the model summary
    xgrid_integration = op.numpy_to_input(
            xgrid_integration, name="integration_grid", custom_shape=(None, 1))

    # 1c Get the original grid
    x_original = Lambda(lambda x: op.op_gather_keep_dims(x, -1, axis=-1), name="x_original_integ")(xgrid_integration)

    # 2. Divide the grid by x depending on the flavour
    x_divided = xDivide()(x_original)

    # 3. Prepare the pdf for integration by dividing by x
    pdf_integrand = Lambda(op.op_multiply, name="pdf_integrand")([x_divided, pdf_xgrid_integration])

    # 4. Integrate the pdf
    pdf_integrated = xIntegrator(weights_array, input_shape=(nx,))(pdf_integrand)

    # 5. Compute the normalization factor
    normalization_factor = MSR_Normalization(output_dim, mode, name="msr_weights")(pdf_integrated)

    # 6. Apply the normalization factor to the pdf
    pdf_normalized = Lambda(lambda pdf_norm: pdf_norm[0] * pdf_norm[1], name="pdf_normalized")([pdf_x, normalization_factor])

    inputs = {
        "pdf_x": pdf_x,
        "pdf_xgrid_integration": pdf_xgrid_integration,
        "xgrid_integration": xgrid_integration,
        }
    model = MetaModel(inputs, pdf_normalized, name="impose_msr")

    return model, xgrid_integration

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

