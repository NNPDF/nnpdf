# n3fit Architecture

## Self-update rule
After any prompt that required reading 2 or more files, append new findings to this file.
Update existing sections in place rather than duplicating. Keep entries concise.

## Key files
- `model_gen.py`: PDF model generation. `_pdfNN_layer_generator` builds the full model; `_generate_nn` builds one replica's NN. `ReplicaSettings` dataclass holds per-replica config including BNN params.
- `backends/keras_backend/base_layers.py`: Layer definitions. `VBDense` = mean-field Bayesian dense (LRT training, weight-sample eval). `layers` dict maps string names to classes.
- `layers/losses.py`: `LossInvcovmat` (chi2), `LossKL` (BNN KL term, reads vb_layers directly), `LossRepulsion` (SVGD repulsion).
- `bnn_wrapper.py`: canonical `BNNPredictor`, `get_vb_layers`, `set_model_eval` for the non-Laplace path. `get_vb_layers` matches `VB_LAYER_CLASSES = (VBDense, CorrelatedLowRankVBDense)`, so it also picks up the correlated layer. `generate_bnn_replica` transfers the posterior by zipping `trainable_weights` and assigning by `.name` (NOT by calling `copy_vb_posterior`, which is defined in this file but currently unused/dead).
- `bnn_wrapper_copy.py`: extended `BNNPredictor` with linearized-Laplace support. Adds `eval_map_model`, `FrozenSampleModel`, `compute_jacobian`, `compute_sigma_theta`, `compute_sigma_f`, `sample_function_space`, `linearized_laplace_samples`. Runcard key `sampling_space: weight|function` selects the path. Also supports `VB_LAYER_CLASSES = (VBDense, CorrelatedLowRankVBDense)` (ported from `bnn_wrapper.py` on 2026-09-18). Here `copy_vb_posterior` (used by both `eval_map_model` and `generate_bnn_replica_from_weights`, unlike in `bnn_wrapper.py`) excludes `random`/`random_b`/`getattr(_, "random_eps", None)` — `random_eps` is `CorrelatedLowRankVBDense`'s extra low-rank noise buffer; omitting it from the skip-list would collapse every pseudo-replica's low-rank component onto the parent's frozen draw, the same class of bug the `random`/`random_b` exclusion already guards against.
- `model_trainer.py`: `ModelTrainer` drives fitting; imports `LossKL`, `LossRepulsion`, `get_vb_layers`.
- `performfit.py` / `performfit_copy.py`: both define a function named `performfit` (same reportengine action name — `n3fit_exec.py::N3FIT_PROVIDERS` picks which module provides it). Both compute `is_bnn = any(layer in {'VBDense', 'VBDense_correlated'} for layer in layer_type)` (kept in sync as of 2026-09-18). `performfit_copy.py`'s BNN path uses `BNNPredictor` from `bnn_wrapper_copy.py` and 0-indexes pseudo-replicas (SLURM-parallel scheme, see below); `performfit.py`'s BNN path uses `BNNPredictor` from `bnn_wrapper.py` and writes 1-indexed replicas.
- `scripts/n3fit_exec.py`: CLI entry point. Parses replica number, builds `N3FitEnvironment`, calls `performfit` via validphys/reportengine. `N3FIT_PROVIDERS`'s first entry (`n3fit.performfit` vs `n3fit.performfit_copy`) has been flipped back and forth locally (uncommitted) — always check this line directly rather than trusting a memory of which one is "active"; don't assume.
- `io/writer.py`: `WriterWrapper` + `storefit`. `WriterWrapper.write_data(replica_path, fitname, save)` iterates replica indices and writes `.exportgrid`, `.json`, `chi2exps.log`. `storefit` converts EVOL-basis model output → LHA-basis YAML exportgrid via `evln2lha`. `XGRID` constant defined here.
- `vpinterface.py`: `N3PDF(pdf_models, fit_basis, Q)` — validphys-compatible wrapper around MetaModel(s). `__call__(xarr, flavours="n3fit")` returns raw EVOL-basis output. Used by `storefit`.
- `mc_dropout_inference.py`: standalone post-training script for MC-dropout fits. Pattern to follow for `bnn_inference.py`.
- `bnn_inference.py`: standalone post-training BNN resampling script. Loads `weights.weights.h5` from an existing pseudo-replica, rebuilds the architecture from the runcard, and draws new weight-space samples without retraining.

## Model tensor shapes
- Input: `(1, n_x, 2)` — batch=1, x-points, [x, log(x)]
- After NN hidden layers: `(1, n_x, n_hidden)`
- NN output (per replica, pre-preprocessing): `(1, 1, n_x, n_flavours)` with replica axis added
- After stacking all replicas: `(1, n_replicas, n_x, n_flavours)`
- Final PDF: `(1, n_replicas, n_x, 14)` in FK basis

## BNN mechanics
- `VBDense` stores `mu_w`, `logsig2_w`, `bias`; training uses LRT (sample activations), eval samples weights once and freezes
- KL is analytic: `KL[N(mu, sigma^2) || N(0, 1/prior_prec)]`
- `LossKL` calls `layer.kl_loss()` directly (stateful pattern — reuse this for new losses)
- `kl_beta` is a `tf.Variable` on `MetaModel`, annealed during training
- `set_model_eval(model)` calls `layer.eval()` on all `VBDense` layers; does NOT delete `mu_w`/`logsig2_w` — they survive `save_weights`.
- Pre-existing bug fixed: `Kops.log(python_float)` (keras.ops) ignores `keras.config.floatx()` and always returns float32, which crashed `_gaussian_kl`'s subtraction against a float64 `logsig2` the moment `kl_loss()` was ever actually called end-to-end. Use `math.log` for plain-Python-float constants combined with float64 tensors, not `Kops.log`.
- Pre-existing bug fixed: `ReplicaSettings` had no `std_init` field even though `model_trainer.py` always passed `std_init=...` into its constructor — `hyperparametrizable` could not build *any* BNN model (flow or not) before this was added.
- `n3fit_exec.py::N3FIT_PROVIDERS` selects which controller runs the fit: `n3fit.performfit` (legacy) vs `n3fit.performfit_copy` (SLURM-parallel BNN controller, 0-indexed pseudo-replica numbering: BNN model `k` writes replicas `k*n_bnn_samples .. (k+1)*n_bnn_samples-1`). Whichever is listed there is what actually runs — check this file first if BNN behavior seems to not match what you read in `performfit.py`.
- `BNNPredictor(pdf_model, n_bnn_samples=3, sampler="weight")` (`bnn_wrapper_copy.py`) has two modes: `sampler="weight"` → `generate_bnn_replica_from_weights()` (classical, self-contained, needs nothing extra); `sampler="function"` → `linearized_laplace_samples()` (linearized Laplace / function-space, needs `pdf_model.invcovmat_per_dataset`/`xgrid_per_dataset`/`obs_wrappers_per_dataset`, all attached by `model_trainer.py::hyperparametrizable` right after `pdf_model` is built, sourced from `_generate_replica_losses`'s per-dataset loop). Runcard key: `parameters.sampling_space: weight|function` (default `weight`).
- Fixed: the Laplace sampler's `_compute_prediction_jacobian_d` feeds `map_model`'s PDF output through `ObservableWrapper._generate_experimental_layer` (FK convolution) to get training-space Jacobians — but `map_model` (built via `single_replica_generator`, `replica_axis=False`) outputs rank-3 `(batch, xgrid, flavours)`, while `Observable.call` (DIS/DY) requires rank-4 `(batch, replicas, xgrid, flavours)`. Fixed by reinserting a size-1 replica axis (`pdf_out[:, tf.newaxis, :, :]`) before the FK convolution.
- `sample_function_space`'s progressive-jitter Cholesky loop (10⁻⁶ → 10⁻¹) logs several benign `Eigen::LLT failed` TF warnings before it finds a jitter that works — expected, not a crash, given `Basic_runcard_qed`'s ill-conditioned `G` (condition number ~10¹¹ even after prior-precision stabilisation, on this toy single-dataset setup).

## NF-enriched VBDense posterior (flow on the last Bayesian layer)
- `backends/keras_backend/normalizing_flow.py`: `RealNVPCouplingLayer`/`RealNVPFlow`, ported from `myScripts/bayes_nf_integral.py`. Identity-initialised (zero-init `s2`/`t2` kernels), per-sample log-det. Do NOT import these into `base_layers.py` directly — that module shadows `Dense` with its own `MetaLayer`-mixed subclass.
- `VBDense(use_flow=True, n_flows=3, flow_hidden=24)`: opt-in per instance, off by default (byte-for-byte unchanged when `use_flow=False`). Flow acts on the WHOLE flattened `(out,in)` weight matrix (`dim=out*in`), not per-row.
- LRT (`_forward_sample_activations`) is invalid once `q(w)` is flow-transformed (it needs Gaussian pre-activations). Flow-enabled layers use `_forward_sample_weights_stochastic` at training time instead: explicit reparam, fresh noise every call, real weight sample `w=flow(z)` materialized and used directly in the matmul.
- `kl_loss()` branches: non-flow layers keep the analytic `_gaussian_kl`; flow layers use `_mc_kl_flow()`, a single-sample MC estimate `log q(w) - log p(w)` computed from `self._last_z/_last_w/_last_logdet`, cached during the most recent stochastic forward pass. This caching is a stateful side-channel: `LossKL` is wired in `model_trainer.py::_pdf_injection` as `f(pdf_layers[0])` with the argument value ignored — the cache is only fresh because that call is topologically downstream of the VBDense forward pass in the same functional-graph trace. Do not reorder that wiring without preserving the dependency.
- Eval/pseudo-replica path (`_forward_sample_weights`, frozen `self.random`) also pushes through the flow, so post-training sampling uses the trained-Gaussian-plus-flow posterior automatically.
- `bnn_wrapper.py`/`bnn_wrapper_copy.py`: `copy_vb_posterior(parent_vb, child_vb)` replaces the old 3-line hardcoded `mu_w`/`logsig2_w`/`bias`-only copy (which silently dropped `bias_logsig2` and would've dropped flow weights too). Generic positional copy over `.weights`, explicitly skipping `random`/`random_b` so each child replica keeps its own independent frozen noise draw.
- Runcard fields: `bayesian_flow` (per-layer list, same convention as `bayesian_bias`), `flow_n_couplings`, `flow_hidden_units` (scalars). Threaded through `ReplicaSettings` → `_generate_nn` → `bnn_inference.py`'s architecture reconstruction.
- Verified (see scratchpad `verify_vbdense_flow.py`/`verify_smoke_fit.py` from the implementing session): gradients reach every flow variable except `s1`/`t1` at step 0 exactly (expected — `s2`/`t2`'s zero-init kernel *is* the backward Jacobian through them, resolves itself from step 1 onward); MC-KL matches analytic KL within MC error at flow=identity; full `_pdfNN_layer_generator` build + training loop with the runcard's real `['dense','VBDense','VBDense']` architecture runs clean (no NaNs, loss decreases, `bayesian_flow=False` reproduces the pre-flow architecture).

## SLURM-parallel BNN training (`performfit_copy.py`)
- `n3fit <runcard> $SLURM_ARRAY_TASK_ID` where task IDs run 1..n_bnn_models
- SLURM task R → `bnn_idx = R - 1` (0-indexed BNN model)
- BNN model `k` writes pseudo-replicas `k*n_bnn_samples .. (k+1)*n_bnn_samples - 1`
- Seed variation: `nnseed + bnn_idx` per model, ensuring independent posteriors
- Bounds check: raises `ValueError` if task ID > `n_bnn_models`
- Non-BNN path unchanged: task R → `replica_R` output

## BNN resampling (`bnn_inference.py`)
- Prerequisite: n3fit must have been run with `save: weights.weights.h5` in the runcard
- All pseudo-replicas from the same BNN model share `mu_w`/`logsig2_w`; loading any one of them is sufficient to resample
- Weight-discovery: scans `replica_{bnn_idx*n_bnn_samples_orig + k}/weights.weights.h5` for k in range(n_bnn_samples_orig)
- Writes `.exportgrid` via `storefit` (needs `validphys.loader.Loader().check_theoryID(theoryid)`)
- Copies source `.json` and `chi2exps.log` to new replica dirs for postfit compatibility
- `--bnn-replica N` (1-indexed) resamples only model N; omit to resample all models
- Both sampling modes supported post-hoc: `--sampler weight` (self-contained) and `--sampler function`
  (linearized Laplace). `function` mode calls `_reconstruct_training_data`, which resolves
  `experiments_data`/`replicas_nnseed_fitting_data_dict`/`posdatasets_fitting_pos_dict`/
  `integdatasets_fitting_integ_dict`/`basis`/`fitbasis` directly via
  `reportengine.resourcebuilder.ResourceBuilder` + `templateparser.string_to_target("datacuts::theory::fitting <name>")`
  (namespace string matches `n3fit_exec.py`'s `FIT_NAMESPACE`) using `N3FitConfig`/`N3FitEnvironment`,
  WITHOUT running the `performfit` action or `N3FitApp.run()` (no `init_output()` side effects — a
  handful of `Environment` attrs it needs, e.g. `replicas`/`output_path`/`replica_path`/`hyperopt`,
  are set by hand instead). The resolved data feeds a real (untrained) `ModelTrainer` — only
  `_generate_static_data`+`_generate_replica_losses` run, no NN is built/trained — to get
  `invcovmat_per_dataset`/`xgrid_per_dataset`/`obs_wrappers_per_dataset`, attached onto the
  weights-loaded `pdf_model` before `BNNPredictor(sampler="function")`. Needs `replica_cli = bnn_idx+1`
  (matches the CLI replica number originally used to train that BNN model, since that's what seeds the
  tr/vl mask) — NOT always `replicas=[1]`. Needs the same theory/commondata on disk the original fit
  used; does NOT need `vp-setupfit`'s output directory (`use_cuts='internal'` recomputes cuts fresh).

## Adding new losses
Follow the `LossKL` pattern: store state in the layer, read it in a `MetaLayer` subclass called with a dummy input. Wire it into `_pdf_injection` in `model_trainer.py`.

## `CorrelatedLowRankVBDense` ("VBDense_correlated") — correlated posterior
- `backends/keras_backend/base_layers.py`: `CorrelatedLowRankVBDense` — `q(w) = N(mu, D + U U^T)`, `rank=0` reproduces `VBDense` exactly. Same interface as `VBDense` (`.mu_w`, `.bias`, `.prior_prec`, `.bayesian_bias`, `.random`/`.random_b`, `.enable_map()`/`.disable_map()`/`.eval()`/`.reset_random()`), PLUS an extra non-trainable low-rank noise buffer `.random_eps` (only present when `rank > 0`) that has no `VBDense` counterpart — see next bullet.
- `model_gen.py::_generate_nn`'s `layer_generator` has a `"VBDense_correlated"` branch (sibling to `"VBDense"`) passing `rank`/`u_init` through; `ReplicaSettings`/`_generate_nn` both carry `rank: int = 4`/`u_init: float = 1e-4` fields alongside (not instead of) the unrelated `bayesian_flow`/`flow_n_couplings`/`flow_hidden_units` fields — these two feature sets (flow-enriched posterior vs. correlated-rank posterior) were added independently and landed as a real merge conflict (resolved 2026-09-18 by keeping both, since they're orthogonal). `model_trainer.py`'s `ReplicaSettings(...)` construction reads `rank`/`u_init` from `params.get(...)` the same way.
- Any code that does `isinstance(layer, VBDense)` to detect "is this a Bayesian layer" must use `VB_LAYER_CLASSES = (VBDense, CorrelatedLowRankVBDense)` instead, or it will silently skip correlated layers. Both `bnn_wrapper.py` and `bnn_wrapper_copy.py` do this in `get_vb_layers`.
- **Any code that copies posterior params parent→child replica must explicitly exclude `random_eps`** (in addition to `random`/`random_b`), or the low-rank noise component collapses across every pseudo-replica of a correlated-layer BNN (identical bug class to the documented `random`/`random_b` regression). `bnn_wrapper.py` avoids the issue structurally (copies only `trainable_weights`, which excludes all three non-trainable buffers automatically); `bnn_wrapper_copy.py`'s `copy_vb_posterior` excludes it explicitly via `getattr(parent_vb, "random_eps", None)`.
- `performfit.py`/`performfit_copy.py`'s `is_bnn` check must test membership in `{'VBDense', 'VBDense_correlated'}`, not equality with `'VBDense'`.
- The Laplace/function-space sampler (`bnn_wrapper_copy.py`'s `compute_jacobian`/`compute_sigma_theta`/etc.) assumes an isotropic Gaussian prior (`prior_prec * I`) and was NOT extended for the correlated layer's `U U^T` prior structure — this is a known gap, not yet hit by any runcard in this repo (`Basic_runcard_bayes.yml` uses plain `VBDense`), left as-is when porting correlated-layer support elsewhere.

## Testing setup (BNN smoke fits)
- The `n3fit` CLI installed in `environment_nnpdf` resolves to a **separate git checkout** (`.../site-packages/n3fit`, remote `ramonpeter/BayesianPDF.git`) — NOT this dev repo. Dev-repo edits have no effect on plain `n3fit` runs. To actually test dev-repo code: `PYTHONPATH=/home/daksh/nnpdf/nnpdf/n3fit/src n3fit <runcard> <replica>` (verify with `python -c "import n3fit; print(n3fit.__file__)"`).
- Quick correctness check for BNN logic changes: copy `runcards/examples/Basic_runcard_bayes.yml` to a scratch dir, drop `parameters.epochs` to ~20, set `stopping_patience: 1.0`, add top-level `parallel_models: false` (required whenever `layer_type` isn't all `'dense'`, else `checks.py::check_consistent_parallel` fails fast). Optionally set `parameters.n_bnn_models`/`n_bnn_samples`/`sampling_space` explicitly (defaults 1/3/`weight`).
- Verified 2026-09-18: both the `performfit.py`+`bnn_wrapper.py` path and the `performfit_copy.py`+`bnn_wrapper_copy.py` path (after porting `VB_LAYER_CLASSES`/`BAYESIAN_LAYER_TYPES` support into the latter, see above) run this smoke config to completion and produce distinct (non-collapsed) per-replica `.exportgrid` values.

## Writer / output file structure
- Per replica: `nnfit/replica_{N}/{fitname}.exportgrid`, `{fitname}.json`, `chi2exps.log`
- Weights (if saved): `nnfit/replica_{N}/weights.weights.h5` — full Keras weights including VBDense variational params
- `WriterWrapper._write_weights` saves `pdf_objects[i]._models[0]` (the inner MetaModel from `N3PDF`)
- `storefit` needs: an `N3PDF`-callable, the replica index, an output path, and a theory object
- Theory object loaded via `validphys.loader.Loader().check_theoryID(theoryid_int)`

## N3PDF / vpinterface
- `N3PDF(pdf_models_or_single, fit_basis=basis, Q=q0)` — wraps one or more MetaModel objects
- `N3PDF.__call__(xarr, flavours="n3fit")` → delegates to `N3LHAPDFSet.__call__` → raw EVOL-basis output
- `_write_weights` accesses `pdf_objects[i]._models[0]` to get the inner MetaModel

## `bnn_wrapper_copy.BNNPredictor` attributes set at training time
- `invcovmat_per_dataset`, `xgrid_per_dataset`, `obs_wrappers_per_dataset` are attached to the MetaModel by `ModelTrainer` during training
- These are only needed by the Laplace sampler; `__init__` now uses `getattr(..., None)` so the class is safe to instantiate from a reloaded model

## `plt_corr.py` (weight-correlation study)
- Standalone analysis script (h5py/numpy/matplotlib only, no TF). Correlates flattened weights across replicas from `<results>/<fit>/nnfit/replica_*/weights.weights.h5`; docs in `~/plt_corr_out/plt_corr_methods.md`. `--last-layer` restricts to the last layer.
- h5 layout: `layers/meta_model/layers/meta_model/layers/{dense,dense_1,..,vb_dense}/vars/<i>`. VBDense vars = bias, mu_w, logsig2_w, [bias_logsig2], random_w, [random_b] (6 vars if bayesian_bias, else 4); frozen eval weight is `mu + sqrt(exp(logsig2))*random_w`. VBDense kernel is `(out,in)`, dense kernel `(in,out)` -> transpose VB to align.
- Samples of one BNN model share mu/logsig2 AND all non-Bayesian layers, so those weights have only `n_bnn_models` distinct values across replicas: correlations of shared weights have effective N = n_models, not N (noise floor ~ 1/sqrt(n_models)). `bay_1x1000` (1 model) has constant hidden layers -> undefined correlation. Results dir names: `bay_AxB` = A models x B samples (bay_10x100 is really 12x100; 1200 replicas each).
- `nnpdf40-like` (fit in results dir) has NO saved weights; the weights fit is `nnpdf40-like_weights` (runcard `runcards/examples/OTHERS/`, jobs in `~/jobs/SNN_JOBS/nnpdf40-like_weights/`). SLURM `MaxArraySize` is 1001, so 1200-replica arrays must be split.
- `plt_corr.py --per-chain` correlates within BNN chains (chain size = `n_bnn_samples` in `<results>/<fit>/filter.yml`; replica r -> chain r // K). `plt_pdf_corr.py` compares PDF-space correlations from each replica's `<fit>.exportgrid` (196 x-points x 14 flavour-basis columns, `xgrid` entries parse as strings in yaml). `bay_100x10` chain 25 (replicas 250-259) is a failed chain (~100 MAD-sigma PDF excursions) -> use `--drop-outliers 50`.
- `plt_corr.py --single-replica R` adds a per-runcard "structure map" z_i*z_j (z = standardised single-replica weight vector; NOT a correlation, deterministic, no ensemble needed) -- useful because the real rho(w,w') across-replica formula is 0/0 undefined for N=1. Caveat: standardising over the whole mixed weight vector lets the largest-scale layer (last VBDense) dominate the colour range, making earlier dense layers look artificially flat. `replica_0`/`replica_1` are NOT central/average replicas in any fit here (BNN pseudo-replicas are 0-indexed with no averaging; standard fits are 1-indexed, no replica_0, no postfit average produced in these runs) -- any replica index is an arbitrary example.
- results/nnpdf40-like_weights was symlinked to the live runcards/examples/OTHERS/nnpdf40-like_weights fit dir (still running, ~30 replicas/batch under SLURM's %30 array throttle) so plt_corr.py/plt_pdf_corr.py can analyse whatever has finished so far without waiting or copying; corr.sh's `cp -r` is now a no-op given the symlink.
