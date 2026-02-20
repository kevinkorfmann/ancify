# Changelog

## 1.5.0 (2026)

Documentation: tutorials for ancestral FASTAs and VCF polarization, plus a species catalogue.

### Documentation

- **Tutorials.** New **Tutorial 4: Getting Ancestral FASTA Files** — how to locate, verify, and read per-chromosome output; sanity checks; loading all chromosomes; confidence encoding recap. New **Tutorial 5: Polarizing VCF Variants** — annotating REF/ALT as ancestral/derived with cyvcf2 and scikit-allel, a complete `polarize_vcf.py` script (AA/DAF INFO fields), unfolded SFS example, and tips on confidence filtering.
- **Species guide.** New **Species catalogue** with 115 commonly studied species: suggested inner and outer outgroups, UCSC assembly identifiers, and approximate divergence times. Grouped by clade (primates, rodents, carnivores, ungulates, birds, fish, Drosophila, plants, fungi, etc.). Includes a short section on generating your own net AXT alignments with lastz when UCSC data is unavailable.

---

## 1.4.0 (2026)

Likelihood-based ancestral reconstruction and expanded installation docs.

### Likelihood method (Felsenstein pruning)

- **New `method: likelihood`.** Infers ancestral alleles using Felsenstein's pruning algorithm on a user-supplied tree with branch lengths. Root posterior probabilities are computed under a continuous-time substitution model and mapped to the same case-encoded confidence scheme (uppercase / lowercase / `n` / `N`).
- **Substitution models.** Four models are supported: **JC69**, **K80**, **HKY85**, and **GTR**. All use normalised rate matrices (one expected substitution per unit branch length). Transition probabilities use `scipy.linalg.expm` (JC69 has a closed-form shortcut).
- **New `ancify.likelihood` module.** `SubstitutionModel` base class, `JC69`, `K80`, `HKY85`, `GTR` classes, `felsenstein_pruning()`, `call_ancestral_base_likelihood()`, `_call_chromosome_likelihood()` worker, and `build_model()` factory.
- **New config fields:** `substitution_model`, `model_kappa`, `model_base_freqs`, `model_rates`, `likelihood_high_threshold`, `likelihood_low_threshold`. Validation requires a tree with leaf names matching outgroups; GTR requires six `model_rates`; base frequencies must sum to ~1.
- **Core dependency:** added **SciPy** (>=1.7) for matrix exponentiation.

### Documentation and README

- **Algorithm docs.** New section on the likelihood method: substitution models, Felsenstein pruning steps, worked example, and comparison with voting/parsimony/ML. Summary section updated for all four methods.
- **Configuration docs.** Field reference and method table updated for `likelihood`; new subsection "Likelihood YAML" with branch-length requirement and GTR example; validation and config recipes updated.
- **Landing page (docs/index.rst).** Likelihood added to the method list and "Why ancify?" bullet.
- **README.** Major expansion: **Installation** now includes prerequisites, core install, optional extras table (`evaluate`, `fast`, `ml`, `docs`, `dev`, `all`), GPU acceleration (voting only), verify-install commands, and quick-reference table. README also updated for four methods throughout (intro, confidence encoding, key fields, How it works, CLI, project structure).

### Tests

- **New `tests/test_likelihood.py`.** Rate-matrix properties (row sums, detailed balance, normalisation), transition-probability properties (P(0)=I, rows sum to 1, equilibrium limit), Felsenstein pruning and posteriors, confidence encoding, agreement with parsimony on unambiguous cases.
- **`tests/test_config.py`.** New `TestValidateLikelihood` for tree requirement, model name, GTR rates, base freqs, and thresholds.
- **`tests/test_ancestral.py`.** New `TestCallAncestralBaseLikelihood` mirroring parsimony tests.

### Version and CLI

- **Version** set to **1.4.0** in `pyproject.toml`.
- **CLI template** (`EXAMPLE_CONFIG` in `cli.py`) updated with commented likelihood example and branch-length tree.

---

## 1.3.0 (2026)

Machine learning-based ancestral calling and documentation updates.

- **ML-based ancestral calling.** New `method: ml` option uses a LightGBM gradient-boosted classifier trained on per-position features (outgroup agreement, GC content, CpG flag, etc.) to predict ancestral alleles. Confidence is derived from predicted class probabilities. Install with `pip install 'ancify[ml]'` (requires `lightgbm` and `scikit-learn`).
- **New `ancify.ml` module.** Feature extraction (`extract_features()`), model loading, and vectorized prediction for full-chromosome runs. Integrates with the existing pipeline via config and CLI.
- **New config field:** `method` now supports `"voting"`, `"parsimony"`, and `"ml"`. For `method: ml`, optional `model_path` points to a trained LightGBM model (or uses a bundled default when available).
- **CLI and config** updated to pass method selection and ML options through to the calling phase.
- **Documentation:** algorithm page and configuration reference updated for the ML method; GPU logo and conf tweaks.
- **Tests:** `tests/test_ml.py` for feature extraction, prediction shape, and integration.
- **Lock file:** `uv.lock` added for reproducible installs.

## 1.2.0 (2026)

Fitch parsimony for tree-based ancestral inference.

- **Fitch parsimony method.** New `method: parsimony` option uses the Fitch (1971) algorithm on a user-supplied Newick phylogenetic tree to infer ancestral alleles. This resolves many positions that the two-tier voting method marks as "unresolved" by leveraging the tree topology.
- **Newick tree parser.** Built-in recursive-descent parser for Newick-format trees (`ancify.parsimony`). Supports branch lengths, quoted labels, and multifurcations.
- **New config fields:** `method` (`"voting"` / `"parsimony"`) and `tree` (inline Newick string or path to `.nwk` file).
- **Config validation** checks that tree leaf names match outgroup species names when parsimony is selected.
- **New `call_ancestral_base_parsimony()` function** in `ancify.ancestral` for programmatic per-position Fitch calls.
- **Comprehensive test suite** for the parsimony module: Newick parsing, Fitch bottom-up/top-down passes, full algorithm with ILS scenarios, missing data handling, and confidence encoding.
- **Documentation updates:** algorithm page with Fitch walkthrough, configuration reference with parsimony examples, API reference for the new module, and updated README.

## 1.1.0 (2026)

GPU acceleration and a vectorized compute backend for much faster Phase 1 and Phase 2 runs.

### Performance

- **GPU-accelerated ancestral calling (Phase 2).** Ancestral state inference runs as a small number of tensor operations on the GPU instead of per-position Python loops. On an NVIDIA A100, the full human genome completes in under 2 minutes (vs. hours on the original scalar path).
- **Vectorized coordinate projection (Phase 1).** Net AXT projection uses NumPy vectorized scatter (CPU) or PyTorch scatter on CUDA. The per-character Python loop is removed, giving roughly 20–50× speedup on CPU.
- **Multi-GPU support.** When using the GPU backend, chromosomes are distributed round-robin across available NVIDIA GPUs. Use the `gpu_devices` config field to restrict which devices are used.
- **Faster gzip decompression.** Optional `isal` (Intel ISA-L) dependency provides 2–5× faster gzip decompression for large AXT files. Install with `pip install 'ancify[fast]'`.

### New module and config

- **`ancify.backend` module.** Central abstraction for CPU/GPU execution: `detect_backend()`, `get_available_gpus()`, `open_gz()`, and vectorized implementations of majority vote, ancestral calling, and block scatter for projection.
- **New config fields:** `backend` (`"auto"` / `"cpu"` / `"gpu"`) and `gpu_devices` (optional list of GPU IDs, e.g. `[0, 1, 2]`). With `backend: auto`, ancify uses the GPU when PyTorch and CUDA are available, otherwise the vectorized CPU path.

### Correctness and docs

- **Bit-identical output.** Vectorized and GPU code paths produce the same results as the original scalar implementation. Tie-breaking, `min_inner_freq` / `min_outer_freq` behaviour, and case-encoded confidence are unchanged.
- **New documentation page:** {doc}`performance` with GPU setup, supported hardware, architecture overview, and tuning tips.

## 1.0.0 (2026)

Initial release.

- Config-driven YAML pipeline for any focal species.
- Three-phase workflow: project, call, evaluate.
- Two-tier inner/outer outgroup voting with case-encoded confidence.
- Support for arbitrary numbers of inner and outer outgroup species.
- Parallel execution via `ProcessPoolExecutor`.
- Optional evaluation against reference ancestral sequences and VCF data.
- CLI with subcommands: `init`, `project`, `call`, `evaluate`, `run`.
- Python API for programmatic use.
- 108 unit and integration tests.
- Example configs for human, mouse, Drosophila, and Brassica rapa.
- Comprehensive documentation with population genetics background, tutorials, and algorithm deep dives.
- Installable with pip or uv.
