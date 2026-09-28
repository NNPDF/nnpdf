"""
Wrapper for BNN (Bayesian Neural Network) inference

This module provides utilities for:
1. Detecting if a model is a BNN (has VBDense layers)
2. Generating pseudo-replicas from BNN weight samples for plotting/analysis
3. Implementing linearized Laplace approximation to sample in function space directly

Linearized Laplace pipeline
----------------------------
There are two distinct Jacobians that must not be confused:

  J_train : (N_data x |theta|)
      Jacobian at experimental kinematic points X_train.
      Used ONLY to build the Gauss-Newton matrix and invert it:
          G         = J_train^T C_exp^{-1} J_train + prior_prec * I
          Sigma_theta = G^{-1}                          [weight space]

  J_pred : (N_grid x |theta|)
      Jacobian at the fine x-grid where PDF bands are wanted.
      Used ONLY to push Sigma_theta forward to function space:
          Sigma_f = J_pred Sigma_theta J_pred^T         [function space]

Sampling is then done directly from N(f_MAP, Sigma_f) — weights are
never sampled, so flat directions (fibers) are automatically annihilated
by the Jacobian projection.

Schematically:

  X_train --> compute_jacobian --> J_train
                                      |
                                compute_sigma_theta <-- C_exp, prior_prec
                                      |
                                  Sigma_theta   (weight space)
                                      |
  X_pred  --> compute_jacobian --> J_pred
                                      |
                                compute_sigma_f
                                      |
                                   Sigma_f      (function space)
                                      |
                               sample_function_space
                                      |
                              f_samples (N_grid x n_samples)
"""

from n3fit.backends.keras_backend.base_layers import CorrelatedLowRankVBDense, VBDense
from n3fit.layers.preprocessing import BayesianPreprocessing
from n3fit.io.writer import XGRID
from n3fit.backends import MetaModel
import tensorflow as tf
import numpy as np

VB_LAYER_CLASSES = (VBDense, CorrelatedLowRankVBDense)


# ---------------------------------------------------------------------------
# Model introspection helpers
# ---------------------------------------------------------------------------

def is_bayesian_model(pdf_model):
    """
    Check if the given pdf_model is a BNN (contains VBDense layers).

    Parameters
    ----------
    pdf_model : MetaModel
        The PDF model to check

    Returns
    -------
    bool
        True if the model is a BNN (has VBDense layers), False otherwise
    """
    return len(get_vb_layers(pdf_model)) > 0


def _get_all_layers_recursively(container):
    """Recursively collect all layers at any depth."""
    layers = [container]
    if hasattr(container, 'layers'):
        for sub_layer in container.layers:
            layers.extend(_get_all_layers_recursively(sub_layer))
    return layers


def get_vb_layers(pdf_model):
    """Return all VBDense layers found anywhere in the model hierarchy."""
    return [
        layer
        for layer in _get_all_layers_recursively(pdf_model)
        if isinstance(layer, VB_LAYER_CLASSES)
    ]


def get_bayesian_preprocessing(pdf_model):
    """Return the BayesianPreprocessing layer, or None if absent."""
    for layer in _get_all_layers_recursively(pdf_model):
        if isinstance(layer, BayesianPreprocessing):
            return layer
    return None


def set_model_eval(replica_model):
    """Set all VBDense and BayesianPreprocessing layers to inference mode."""
    for layer in get_vb_layers(replica_model):
        layer.eval()
    preproc = get_bayesian_preprocessing(replica_model)
    if preproc is not None:
        preproc.eval()


def copy_vb_posterior(parent_vb, child_vb):
    """
    Copy every trained VBDense posterior weight from parent to child, generically
    (positional match over `.weights`, which is deterministic since both layers were
    built by the same code path), EXCLUDING the frozen per-replica noise buffers
    (`random`/`random_b`, and `random_eps` for `CorrelatedLowRankVBDense`) so each
    child replica keeps its own independently-drawn posterior sample instead of
    duplicating the parent's.
    """
    for parent_w, child_w in zip(parent_vb.weights, child_vb.weights):
        if (
            parent_w is parent_vb.random
            or (parent_vb.bayesian_bias and parent_w is parent_vb.random_b)
            or parent_w is getattr(parent_vb, "random_eps", None)
        ):
            continue
        child_w.assign(parent_w)

def eval_map_model(pdf_model):
    """Evaluate the pdf model at MAP point estimate"""
    vb_layers = get_vb_layers(pdf_model)

    # Enable map mode
    for layer in vb_layers:
        layer.enable_map()

    map_model = pdf_model.single_replica_generator(0)
    # Transfer preprocessing alpha / beta. MUST run BEFORE copy_vb_posterior/reset_random
    map_model.set_replica_weights(
        pdf_model.get_replica_weights(0), i_replica=0
    )

    # Transfer variational posterior parameters
    new_vb_layers = get_vb_layers(map_model)
    for parent_vb, child_vb in zip(vb_layers, new_vb_layers):
        copy_vb_posterior(parent_vb, child_vb)
        child_vb.reset_random()

    set_model_eval(map_model)

    # Disable map mode 
    for layer in vb_layers:
        layer.disable_map()

    return map_model


# ---------------------------------------------------------------------------
# BNNPredictor
# ---------------------------------------------------------------------------

class FrozenSampleModel:
    def __init__(self, map_model, f_sample, n_flavours=14):
        self._map_model = map_model
        self._f_sample = tf.constant(f_sample, dtype=tf.float32)
        self._n_flavours = n_flavours
        self._n_x = len(f_sample) // n_flavours

    def __call__(self, inputs=None, training=False):
        return self._f_sample

    def predict(self, x=None, **kwargs):
        # Delegate to map_model for arbitrary x queries (arc length, integrability etc.)
        # Only override when queried at exactly the stored grid size
        result = self._map_model.predict(x, **kwargs)
        if result.shape[1] == self._n_x:
            return self._f_sample.numpy().reshape(1, self._n_x, self._n_flavours)
        return result

    def __getattr__(self, name):
        return getattr(self._map_model, name)

class BNNPredictor:
    """
    Predictor class for BNNs.

    Handles two inference modes:

    1. Weight-space sampling (``generate_bnn_replica_from_weights``)
       Classical approach: sample theta ~ q(theta), run nonlinear forward pass.
       Suffers from flat-direction Sdegeneracy in overparameterised nets.

    2. Linearized Laplace in function space 
       Steps:
         a. ``compute_jacobian(X_train)``       -> J_train
         b. ``compute_sigma_theta(J_train, ...)``-> Sigma_theta  [weight space]
         c. ``compute_jacobian(X_pred)``         -> J_pred
         d. ``compute_sigma_f(J_pred, ...)``     -> Sigma_f      [function space]
         e. ``sample_function_space(X_pred, Sigma_f)`` -> PDF replicas

       Flat directions are annihilated by the Jacobian projection and never
       appear in Sigma_f, regardless of prior width.
    """

    def __init__(self, pdf_model, n_bnn_samples=3, sampler="weight"):
        """
        Parameters
        ----------
        pdf_model : MetaModel
            Trained PDF model with VBDense layers.
        map_model : MetaModel
            PDF model evaluated at MAP point estimates
        n_samples : int
            Default number of samples / replicas to generate.
        """
        self.pdf_model = pdf_model
        self.map_model = eval_map_model(pdf_model)
        self.n_samples = n_bnn_samples
        # These are set on the model during training; use getattr so the class
        # also works when instantiated from a reloaded weights file (bnn_inference.py).
        # They are only accessed by the Laplace (function-space) sampler.
        self.invcovmat_per_dataset = getattr(pdf_model, "invcovmat_per_dataset", None)
        self.xgrid_per_dataset = getattr(pdf_model, "xgrid_per_dataset", None)
        self.obs_wrappers = getattr(pdf_model, "obs_wrappers_per_dataset", None)
        self.vb_layers = get_vb_layers(pdf_model)
        self.prior_prec = self.vb_layers[0].prior_prec
        self.bayesian_preproc = get_bayesian_preprocessing(pdf_model)
        self.sampler = sampler

        # Only mu_w and bias — logsig2_w controls posterior width but not the
        # MAP forward pass, so it is excluded from the Jacobian computation.
        self.parameters = []
        map_vb_layers = get_vb_layers(self.map_model)
        for layer in map_vb_layers:
            self.parameters.append(layer.mu_w)
            self.parameters.append(layer.bias)
        self._param_values = self.parameters

    def pdf_sampler(self):
        if self.sampler == "weight":
            return self.generate_bnn_replica_from_weights()
        elif self.sampler == "function":
            return self.linearized_laplace_samples()
        else:
            raise ValueError(f'Unknown sample space {self.sampler}. Please choose either "weight" or "function"')
    # ------------------------------------------------------------------
    # Mode 1: weight-space replica generation (classical BNN)
    # ------------------------------------------------------------------

    def generate_bnn_replica_from_weights(self):
        """
        Generate replicas by sampling weights from the variational posterior
        and running a full nonlinear forward pass for each sample.

        Returns
        -------
        list of MetaModel
            One model per replica, each with fixed (sampled) weights.
        """
        replica_models = []

        for _ in range(self.n_samples):
            replica = self.pdf_model.single_replica_generator(0)

            # Transfer preprocessing alpha / beta (and, redundantly but harmlessly,
            # the NN's dense/VBDense weights -- copy_vb_posterior below re-copies the
            # VBDense posterior params). MUST run BEFORE copy_vb_posterior/reset_random:
            # set_replica_weights blindly overwrites every weight of the "NN" layer,
            # including VBDense's non-trainable `random`/`random_b` buffers, so doing
            # it after reset_random() clobbers the fresh per-replica noise draw back
            # to the parent's frozen value, collapsing every pseudo-replica of a
            # given BNN model onto the same point regardless of n_samples (this was a
            # real regression: histograms showed one spike per trained BNN model
            # instead of a smooth per-model spread).
            replica.set_replica_weights(
                self.pdf_model.get_replica_weights(0), i_replica=0
            )

            # Transfer variational posterior parameters
            new_vb_layers = get_vb_layers(replica)
            for parent_vb, child_vb in zip(self.vb_layers, new_vb_layers):
                copy_vb_posterior(parent_vb, child_vb)
                # Must be the LAST thing touching this replica's weights (see note above).
                child_vb.reset_random()

            # Fix weights and preprocessing for this replica
            set_model_eval(replica)
            replica_models.append(replica)

        return replica_models

    # ------------------------------------------------------------------
    # Mode 2: linearized Laplace in function space
    # ------------------------------------------------------------------
    # To-do: writer.write_data needs to accomodate direct calculations using pdfs; or 
    # return MetaModel pdfs instead of tf ones --> we do the later with FrozenSampleModel
    # ------------------------------------------------------------------
    def compute_jacobian(self, X):
        """
        Jacobian of MAP model PDF output w.r.t. MAP parameters at raw x-grid X.

        Used to push Sigma_theta to function space: Sigma_f = J Sigma_theta J^T.
        This gives posterior variance over PDF values, NOT over theory predictions.

        Parameters
        ----------
        X : tf.Tensor or np.array, shape (N,)
            x-points to evaluate at.

        Returns
        -------
        J : tf.Tensor, shape (N * n_flavours, n_params), dtype float64
        """
        X_tensor = tf.cast(tf.reshape(X, [1, -1, 1]), tf.float64)

        full_input = self.map_model._parse_input({"pdf_input": X_tensor})
        full_input = {
            k: tf.cast(tf.constant(v), tf.float64) if not isinstance(v, tf.Tensor)
            else tf.cast(v, tf.float64)
            for k, v in full_input.items()
        }

        with tf.GradientTape() as tape:
            f = self.map_model(full_input, training=False)
            f = tf.reshape(tf.cast(f, tf.float64), [-1])  # (N * n_flavours,)

        grads = tape.jacobian(f, self._param_values)
        J = tf.concat(
            [tf.reshape(g, [tf.shape(f)[0], -1]) for g in grads],
            axis=1
        )
        return tf.cast(J, tf.float64)  # (N * n_flavours, n_params)

    def _compute_prediction_jacobian_d(self, x_input_d, obs_wrapper_d):
        """
        Jacobian of training theory predictions for dataset d w.r.t. MAP parameters.

        Differentiates through the full pipeline:
            x_input_d -> MAP model -> PDF -> FK convolution (DIS/DY) -> T_d (masked to training set)

        J[i, j] = dT_i / dtheta_j  where T_i is the i-th training-masked prediction.

        Parameters
        ----------
        x_input_d : np.ndarray, shape (1, N_x_d)
            x-grid for this dataset (from pdf_model.xgrid_per_dataset).
        obs_wrapper_d : ObservableWrapper
            Training observable wrapper for this dataset (from pdf_model.obs_wrappers_per_dataset).
            Its _generate_experimental_layer applies the FK convolution and training mask.

        Returns
        -------
        J : tf.Tensor, shape (N_d_tr, n_params), dtype float64
            N_d_tr = number of training data points after mask.
        """
        x_tensor = tf.cast(tf.reshape(x_input_d, [1, -1, 1]), tf.float64)

        full_input = self.map_model._parse_input({"pdf_input": x_tensor})
        full_input = {
            k: tf.cast(tf.constant(v), tf.float64) if not isinstance(v, tf.Tensor)
            else tf.cast(v, tf.float64)
            for k, v in full_input.items()
        }

        with tf.GradientTape() as tape:
            pdf_out = self.map_model(full_input, training=False)
            # map_model is built with replica_axis=False (single_replica_generator), so its
            # PDF output is rank-3 (batch, xgrid, flavours). Observable.call (DIS/DY, invoked
            # inside _generate_experimental_layer) requires rank-4 (batch, replicas, xgrid,
            # flavours) -- reinsert the size-1 replica axis that was stripped.
            pdf_out = pdf_out[:, tf.newaxis, :, :]
            # _generate_experimental_layer: FK convolution + optional rotation + training mask
            T_d = obs_wrapper_d._generate_experimental_layer(pdf_out)
            T_d = tf.reshape(tf.cast(T_d, tf.float64), [-1])  # (N_d_tr,)

        grads = tape.jacobian(T_d, self._param_values)
        J = tf.concat(
            [tf.reshape(g, [tf.shape(T_d)[0], -1]) for g in grads],
            axis=1
        )
        return tf.cast(J, tf.float64)  # (N_d_tr, n_params)

    def compute_sigma_theta(self):
        """
        Gauss-Newton posterior covariance in weight space.

        Correct formulation for PDF fitting:

            G_data = Sigma_d  J_pred_d^T  C_d^{-1}  J_pred_d

        where J_pred_d[i, j] = dT_i^d / dtheta_j is the Jacobian of training
        theory predictions through the FK table convolution, NOT the Jacobian
        of raw PDF values.  The two differ by FK @ (dPDF/dtheta) vs (dPDF/dtheta)
        
        Eigendecomposing G_data only and adding prior_prec analytically bounds
        all eigenvalues of G from below by prior_prec, so inversion is always
        numerically safe (Immer et al. 2410.16901, VIKING 2510.23684).
        """
        if self.obs_wrappers is None:
            raise RuntimeError(
                "pdf_model.obs_wrappers_per_dataset not set. "
                "Ensure model_trainer stores it after model generation."
            )

        n_params = int(sum(tf.size(p) for p in self._param_values))
        G_data = tf.zeros((n_params, n_params), dtype=tf.float64)

        for i, (invcovmat_d, x_input_d, obs_wrapper_d) in enumerate(
            zip(self.invcovmat_per_dataset, self.xgrid_per_dataset, self.obs_wrappers)
        ):
            C_inv = tf.cast(invcovmat_d[0], tf.float64)           # (N_d_tr, N_d_tr)
            C_inv = 0.5 * (C_inv + tf.transpose(C_inv))           # symmetrise

            # J_d: Jacobian of training predictions through FK table, (N_d_tr, n_params)
            J_d = self._compute_prediction_jacobian_d(x_input_d, obs_wrapper_d)
            print(f"[SAMPLING]: Dataset {i}: J_d shape = {J_d.shape}, C_inv shape = {C_inv.shape}")

            G_data = G_data + tf.transpose(J_d) @ C_inv @ J_d

        G_data = 0.5 * (G_data + tf.transpose(G_data))

        eigenvalues_data, eigenvectors = tf.linalg.eigh(G_data)   # ascending, (n_params,)

        # G eigenvalues = lambda_data + prior_prec >= prior_prec > 0 — safe to invert
        prior_prec = tf.cast(self.prior_prec, tf.float64)
        G_eigenvalues = eigenvalues_data + prior_prec

        cond = (tf.reduce_max(G_eigenvalues) / tf.reduce_min(G_eigenvalues)).numpy()
        print(f"[SAMPLING]: G condition number after prior stabilisation: {cond:.3e}")

        # Null-space dirs get 1/prior_prec (prior uncertainty); image-space get 1/(λ_d + prior_prec)
        Sigma_theta = (
            eigenvectors
            @ tf.linalg.diag(1.0 / G_eigenvalues)
            @ tf.transpose(eigenvectors)
        )                                                          # (n_params, n_params)

        G = G_data + prior_prec * tf.eye(n_params, dtype=tf.float64)
        return Sigma_theta, G

    def compute_sigma_f(self, J_pred, Sigma_theta):
        """
        Push Sigma_theta to function space: Sigma_f = J_pred Sigma_theta J_pred^T.

        J_pred is the Jacobian of raw PDF values at the prediction grid (from
        compute_jacobian), so Sigma_f gives posterior variance over PDF values.

        Parameters
        ----------
        J_pred : tf.Tensor, shape (N_grid_outputs, n_params), float64
        Sigma_theta : tf.Tensor, shape (n_params, n_params), float64

        Returns
        -------
        Sigma_f : tf.Tensor, shape (N_grid_outputs, N_grid_outputs), float64
        """
        J = tf.cast(J_pred, tf.float64)
        S = tf.cast(Sigma_theta, tf.float64)
        return J @ S @ tf.transpose(J)
    
    def sample_function_space(self, f_map, Sigma_f, n_samples=None):
        """
        Sample from N(f_MAP, Sigma_f) directly in function space.

        Parameters
        ----------
        f_map : tf.Tensor, shape (N_grid_outputs,)
        Sigma_f : tf.Tensor, shape (N_grid_outputs, N_grid_outputs)
        n_samples : int

        Returns
        -------
        f_samples : tf.Tensor, shape (N_grid_outputs, n_samples)
        f_map : tf.Tensor, shape (N_grid_outputs,)
        """
        if n_samples is None:
            n_samples = self.n_samples

        N = tf.shape(f_map)[0]

        # Symmetrise first
        Sigma_f = 0.5 * (Sigma_f + tf.transpose(Sigma_f))

        # Try progressive jitter
        L = None
        for exp in range(8):
            jitter = 10**(-6 + exp)
            Sigma_stable = Sigma_f + jitter * tf.eye(N, dtype=tf.float32)
            L = tf.linalg.cholesky(Sigma_stable)                            # (N, N)
            if not tf.reduce_any(tf.math.is_nan(L)):
                print(f"Cholesky succeeded with jitter={jitter}")
                break
            L = None

        if L is None:
            print("Warning: Cholesky failed with jitter, falling back to eigendecomposition")
            # Clip negative eigenvalues to zero
            eigenvalues, eigenvectors = tf.linalg.eigh(Sigma_f)
            eigenvalues = tf.maximum(eigenvalues, 1e-8)
            # Reconstruct L = V @ diag(sqrt(lambda))
            L = eigenvectors @ tf.linalg.diag(tf.sqrt(eigenvalues))

        z = tf.random.normal((N, n_samples), dtype=tf.float32)
        f_samples = f_map[:, None] + L @ z                      # (N, n_samples)
        return f_samples

    def linearized_laplace_samples(self, X_pred=None):
        """
        Full linearized Laplace pipeline.

        Steps:
          1. Build G_data = Σ_d J_pred_d^T C_d^{-1} J_pred_d  (prediction Jacobians through FK)
          2. Invert via eigendecomposition to get Sigma_theta  (weight-space covariance)
          3. Compute J_PDF at prediction x-grid                 (PDF Jacobian, not prediction)
          4. Push to function space: Sigma_f = J_PDF Sigma_theta J_PDF^T
          5. Sample from N(f_MAP, Sigma_f) and wrap as FrozenSampleModel
        """
        from n3fit.io.writer import XGRID

        print("\n[SAMPLING]: Laplace SAMPLING starts")

        if X_pred is None:
            X_pred = XGRID

        # Diagnostic: C_inv condition
        for i, invcovmat_d in enumerate(self.invcovmat_per_dataset):
            C_inv = tf.cast(invcovmat_d[0], tf.float64)
            eigvals = tf.linalg.eigvalsh(C_inv)
            print(f"[SAMPLING]: Dataset {i}: C_inv min eigenvalue = {tf.reduce_min(eigvals).numpy():.6e}")

        # Step 1 — weight-space posterior covariance via prediction Jacobians through FK tables
        Sigma_theta, G = self.compute_sigma_theta()

        eigvals_G = tf.linalg.eigvalsh(G)
        print(f"[SAMPLING]: G condition number: {(tf.reduce_max(eigvals_G)/tf.reduce_min(eigvals_G)).numpy():.3e}")
        print(f"[SAMPLING]: G min eigenvalue: {tf.reduce_min(eigvals_G).numpy():.3e}")

        # Step 2 — Jacobian of raw PDF values at prediction grid (to push Sigma_theta to PDF space)
        J_pred = self.compute_jacobian(tf.cast(X_pred, tf.float64))  # (N_grid * 14, n_params)

        # Step 3 — function-space covariance over PDF values
        Sigma_f = tf.cast(self.compute_sigma_f(J_pred, Sigma_theta), tf.float32)

        # Step 4 — evaluate f_MAP at prediction grid
        X_pred_tensor = tf.cast(tf.reshape(X_pred, [1, -1, 1]), tf.float64)
        full_input = self.map_model._parse_input({"pdf_input": X_pred_tensor})
        full_input = {
            k: tf.cast(tf.constant(v), tf.float64) if not isinstance(v, tf.Tensor)
            else tf.cast(v, tf.float64)
            for k, v in full_input.items()
        }
        f_map = tf.reshape(
            tf.cast(self.map_model(full_input, training=False), tf.float32), [-1]
        )                                  Matthias Beck & Sinai Robins, Computing the Continuous Discretely                        # (N_grid * 14,)

        # Step 5 — sample from N(f_MAP, Sigma_f)
        f_samples = self.sample_function_space(f_map, Sigma_f, self.n_samples)

        print("[SAMPLING]: Laplace SAMPLING ends")

        return [FrozenSampleModel(self.map_model, f_samples[:, s]) for s in range(f_samples.shape[1])]
    

