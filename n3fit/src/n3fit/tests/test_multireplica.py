import numpy as np

from n3fit.model_gen import generate_pdf_model


def test_replica_split():
    """Check that multi replica pdf and concatenated single output pdfs agree"""
    num_replicas = 3
    replica_axis = 1
    fake_fl = [
        {"fl": i, "largex": [0, 1], "smallx": [1, 2]}
        for i in ["u", "ubar", "d", "dbar", "c", "g", "s", "sbar"]
    ]
    pdf_model = generate_pdf_model(
        nodes=[8],
        activations=["linear"],
        seed=0,
        flav_info=fake_fl,
        fitbasis="FLAVOUR",
        num_replicas=num_replicas,
    )
    np.random.seed(0)
    fake_input = {
        'pdf_input': np.random.rand(1, 5, 1),
        'integrator_input': np.random.rand(1, 2_000, 1),
    }

    output_full = pdf_model(fake_input)

    pdf_models = pdf_model.split_replicas()
    output_split = [pdf(fake_input) for pdf in pdf_models]
    output_split_stacked = np.stack(output_split, axis=replica_axis)

    np.testing.assert_allclose(output_full, output_split_stacked, rtol=1e-5)
