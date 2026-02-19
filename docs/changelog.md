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
