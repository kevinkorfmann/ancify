# Changelog

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
- Example configs for human, mouse, and Drosophila.
- Installable with pip or uv.
