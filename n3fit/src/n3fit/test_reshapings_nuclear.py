import numpy as np
import tensorflow as tf

from n3fit.backends import Input, Lambda, MetaModel
from n3fit.backends import operations as op
from n3fit.layers import MSR_Normalization, xDivide, xIntegrator
from n3fit.msr import gen_integration_input

nx = 5
num_unique_As = 2
output_dim = 14

integrator_input = tf.constant(list(range(nx)), dtype=tf.float32, shape=(1, nx, 1))
pdf_input_A_unique = tf.constant([1, 10], shape=(1, num_unique_As, 1), dtype=tf.float32)

print(f"Input shapes, x: {integrator_input.shape}, A: {pdf_input_A_unique.shape}")
x_repeated, A_repeated = op.all_combinations(
    integrator_input, pdf_input_A_unique, N_a=nx, N_b=num_unique_As
)
print(f"Output shapes, x: {x_repeated.shape}, A: {A_repeated.shape}")
print(x_repeated[0, :, 0].numpy())
print(A_repeated[0, :, 0].numpy())

pdf_xgrid_integration = tf.repeat(x_repeated + A_repeated, output_dim, axis=2)
print(f"Shape pdf_xgrid_integration: {pdf_xgrid_integration.shape}")
print(pdf_xgrid_integration[0, :, 0].numpy())

pdf_xgrid_integration_reshaped = Lambda(
    lambda x: op.reshape(x, (1, nx, num_unique_As, output_dim)), name="reshape_As"
)(pdf_xgrid_integration)
print(f"Shape pdf_xgrid_integration_reshaped: {pdf_xgrid_integration_reshaped.shape}")
print(pdf_xgrid_integration_reshaped[0, :, :, 0].numpy())

x_divided = tf.constant([1, 2, 3, 4, 5], shape=(1, nx, 1), dtype=tf.float32)
pdf_integrand = Lambda(
    lambda x_pdf: op.op_multiply([op.batchit(x_pdf[0], batch_dimension=2), x_pdf[1]]),
    name="pdf_integrand",
)([x_divided, pdf_xgrid_integration_reshaped])
print(f"Shape pdf_integrand: {pdf_integrand.shape}")
print(pdf_integrand[0, :, :, 0].numpy())

xgrid_integration, weights_array = gen_integration_input(nx)
weights_array = tf.constant([0.1] * nx, shape=(nx, 1), dtype=tf.float32)
print(weights_array)
pdf_integrated = xIntegrator(weights_array, input_shape=(nx,))(pdf_integrand)
print(f"Shape pdf_integrated: {pdf_integrated.shape}")
print(pdf_integrated[0, :, 0].numpy())

photon_integral = tf.constant([3], shape=(1, 1), dtype=tf.float32)
photon_integrals = Lambda(
    lambda x: op.repeat(op.batchit(x, batch_dimension=2), num_unique_As, axis=1),
    name="photon_integral_copied",
)(photon_integral)
print(f"Shape photon_integrals: {photon_integrals.shape}")

normalization_factor = MSR_Normalization('all', num_unique_As=num_unique_As, name="msr_weights")(
    pdf_integrated, photon_integrals
)
print(f"Shape normalization_factor: {normalization_factor.shape}")
print(normalization_factor.numpy())
