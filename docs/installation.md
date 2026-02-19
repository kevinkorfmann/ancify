# Installation

## Quick install

```bash
git clone https://github.com/kevinkorfmann/ancify.git
cd ancify
pip install .
```

That is all you need for the core pipeline (Phases 1 and 2).

## Install options

### With evaluation dependencies

Phase 3 (evaluation against a reference and/or VCF data) requires additional packages:

```bash
pip install '.[evaluate]'
```

### With development dependencies

For running the test suite:

```bash
pip install '.[dev]'
pytest
```

### With GPU acceleration

For large genomes, GPU acceleration can reduce runtime from hours to minutes.
Install PyTorch with CUDA support and the optional fast gzip decompressor:

```bash
# PyTorch with CUDA 12.x (adjust for your driver version)
pip install torch --index-url https://download.pytorch.org/whl/cu128

# Fast gzip decompression (optional, 2-5x faster Phase 1)
pip install '.[fast]'
```

See {doc}`performance` for full details on GPU setup, supported hardware, and
configuration.

### Everything at once

```bash
pip install '.[all]'
```

### Using uv (faster alternative to pip)

[uv](https://github.com/astral-sh/uv) is a fast Python package installer:

```bash
uv pip install .
# or with extras:
uv pip install '.[all]'
```

## Requirements

| Package | Version | Required? | Purpose |
|---------|---------|-----------|---------|
| Python | >= 3.8 | Yes | Runtime |
| PyYAML | >= 5.0 | Yes | YAML config parsing |
| NumPy | >= 1.20 | Yes | Array operations in projection and calling |
| PyTorch | >= 2.0 | No | GPU-accelerated ancestral calling ({doc}`performance`) |
| isal | >= 1.0 | No | 2–5× faster gzip decompression for Phase 1 |
| scikit-allel | >= 1.3 | No | VCF reading (Phase 3 evaluation) |
| matplotlib | >= 3.0 | No | Plotting (Phase 3 evaluation) |
| pytest | >= 7.0 | No | Running the test suite |

## Verify the installation

```bash
ancify --help
```

Expected output:

```text
usage: ancify [-h] [-v] {init,project,call,evaluate,run} ...

Ancestral allele polarization pipeline.

positional arguments:
  {init,project,call,evaluate,run}
    init                Generate a template configuration file
    project             Phase 1: project outgroup alignments
    call                Phase 2: call ancestral alleles
    evaluate            Phase 3: evaluate calls
    run                 Run all phases

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Enable debug logging
```

You can also invoke it as a Python module:

```bash
python -m ancify --help
```

## Platform notes

ancify is pure Python and works on Linux, macOS, and Windows. However:

- **Linux** is recommended for production runs. The pipeline is I/O-intensive and benefits from fast storage (SSD). GPU acceleration with CUDA is fully supported.
- **macOS** works well. If you use Apple Silicon, ensure your Python and NumPy are ARM-native for best performance. GPU acceleration via MPS is not yet supported.
- **Windows** works but has not been extensively tested. WSL2 is recommended for large runs. GPU acceleration works under WSL2 with NVIDIA drivers.

### Memory requirements

Phase 2 loads all projected sequences for a chromosome simultaneously. Memory usage scales as:

```text
memory ≈ num_cpus × (num_inner + num_outer) × chromosome_length
```

For human chr1 (~249 Mb) with 4 outgroups and 24 parallel workers, peak memory can reach ~24 GB. Reduce `num_cpus` if memory is limited.

## Troubleshooting

### `ModuleNotFoundError: No module named 'yaml'`

PyYAML was not installed. Run:

```bash
pip install pyyaml
```

Or reinstall ancify, which will pull it in automatically:

```bash
pip install .
```

### `ModuleNotFoundError: No module named 'allel'`

You are trying to run Phase 3 evaluation without the evaluation extras:

```bash
pip install '.[evaluate]'
```

### `command not found: ancify`

The `ancify` script was not installed on your PATH. Common fixes:

1. **Check your PATH:** `pip show ancify` will show the installation location.
2. **Use the module form:** `python -m ancify` always works if the package is installed.
3. **Virtual environment:** Make sure your virtual environment is activated.

### `SyntaxError` on older Python

ancify requires Python >= 3.8. Check your version:

```bash
python --version
```

If you need an isolated environment:

```bash
conda create -n ancify python=3.10
conda activate ancify
pip install .
```
