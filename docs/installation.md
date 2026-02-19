# Installation

## From source (pip)

```bash
git clone https://github.com/kevinkorfmann/ancify.git
cd ancify
pip install .
```

## From source (uv)

```bash
git clone https://github.com/kevinkorfmann/ancify.git
cd ancify
uv pip install .
```

## With evaluation extras

The evaluation phase (Phase 3) needs `scikit-allel` and `matplotlib`:

```bash
pip install '.[evaluate]'
# or
uv pip install '.[evaluate]'
```

## With development extras

For running the test suite:

```bash
pip install '.[dev]'
pytest
```

## Everything

```bash
pip install '.[all]'
```

## Requirements

| Package | Version | Purpose |
|---------|---------|---------|
| Python | >= 3.8 | Runtime |
| PyYAML | >= 5.0 | Config parsing |
| NumPy | >= 1.20 | Array operations (evaluation) |
| scikit-allel | >= 1.3 | VCF reading (optional, Phase 3) |
| matplotlib | >= 3.0 | Plotting (optional, Phase 3) |
| pytest | >= 7.0 | Testing (optional, development) |

## Verifying the installation

```bash
ancify --help
```

You should see the CLI help text listing all available subcommands.
