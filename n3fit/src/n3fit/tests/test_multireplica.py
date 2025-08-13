import numpy as np

from n3fit.model_gen import ReplicaSettings, generate_pdf_model

EPS = 1e-9
FAKE_FL = [
    {"fl": i, "largex": [0.5, 1.5], "smallx": [1.5, 2.5]}
    for i in ["u", "ubar", "d", "dbar", "c", "g", "s", "sbar"]
]


def test_replica_split():
    """Check that multi replica pdf and concatenated single output pdfs agree"""
    num_replicas = 3
    replica_axis = 1
    rps = num_replicas * [ReplicaSettings(nodes=[8], activations=["linear"], seed=34)]
    pdf_model = generate_pdf_model(rps, flav_info=FAKE_FL, fitbasis="FLAVOUR")
    rng = np.random.default_rng(seed=34)
    pdf_input = np.maximum(rng.random((1, 5, 1)), EPS)
    int_input = np.maximum(rng.random((1, 2_000, 1)), EPS)

    fake_input = {
        'pdf_input': np.sort(pdf_input, axis=1),
        'xgrid_integration': np.sort(int_input, axis=1),
    }

    output_full = pdf_model(fake_input)

    pdf_models = pdf_model.split_replicas()
    output_split = [pdf(fake_input) for pdf in pdf_models]
    output_split_stacked = np.stack(output_split, axis=replica_axis)

    np.testing.assert_allclose(output_full, output_split_stacked, rtol=1e-5)


def test_multimodel(seed=42, xlen=5):
    """Check that we can run different models, with different settings,
    in one go.

    This tests runs 3 replicas with 1, 2, and 3 layers respectively.
    """
    nodes = [20, 10, 8]
    activations = ["tanh", "sigmoid", "linear"]
    init_array = ["glorot_normal", "glorot_uniform", "random_uniform"]

    rps = []
    for i, initialization in enumerate(init_array):
        idx = i + 1
        rps.append(
            ReplicaSettings(
                nodes=nodes[-idx:],
                activations=activations[-idx:],
                seed=seed + idx,
                initializer=initialization,
            )
        )

    rng = np.random.default_rng(seed=seed)
    pdf_input = np.maximum(rng.random((1, xlen, 1)), EPS)
    int_input = np.maximum(rng.random((1, 2000, 1)), EPS)
    fake_input = {
        'pdf_input': np.sort(pdf_input, axis=1),
        'xgrid_integration': np.sort(int_input, axis=1),
    }

    pdf_model = generate_pdf_model(rps, flav_info=FAKE_FL, fitbasis="FLAVOUR")
    output_full = pdf_model(fake_input)
    # Check that the output size is what we expect
    np.testing.assert_array_equal(output_full.shape, (1, len(rps), xlen, 14))
    # And now check that the split model has the right layers
    single_replicas = pdf_model.split_replicas()

    for i, model in enumerate(single_replicas):
        len(model.get_layer("all_NNs").weights) == 2 * (i + 1)
