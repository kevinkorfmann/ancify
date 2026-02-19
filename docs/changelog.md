# Changelog

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

GPU acceleration and vectorized compute backend.

- **GPU-accelerated ancestral calling.** Phase 2 now runs as ~15 tensor operations on the GPU instead of 248 million Python function calls per chromosome. On an NVIDIA A100, the full human genome completes in under 2 minutes.
- **Vectorized coordinate projection.** Phase 1 replaces the per-character Python loop with NumPy vectorized scatter, yielding 20–50× speedup on CPU alone.
- **Multi-GPU support.** Chromosomes are distributed round-robin across available NVIDIA GPUs. Configure with the new `backend` and `gpu_devices` config fields.
- **Faster gzip decompression.** Optional `isal` dependency (Intel ISA-L) provides 2–5× faster gzip decompression for large AXT files. Install with `pip install 'ancify[fast]'`.
- **New `ancify.backend` module.** Central abstraction for CPU/GPU execution with `detect_backend()`, `get_available_gpus()`, `open_gz()`, and vectorized implementations of majority vote, ancestral calling, and block scatter.
- **New config fields:** `backend` (`"auto"` / `"cpu"` / `"gpu"`) and `gpu_devices` (list of GPU IDs).
- **Bit-identical output.** All vectorized and GPU code paths produce identical results to the original scalar implementation. Tie-breaking, frequency thresholds, and confidence encoding are preserved exactly.
- **New documentation page:** {doc}`performance` with full details on GPU setup, supported hardware, architecture, and tuning tips.

## 1.0.0 (2025)

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
