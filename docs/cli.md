# Command-Line Interface

ancify provides a single CLI command with several subcommands. The package can
be invoked either as `ancify` (after installation) or as `python -m ancify`.

## Synopsis

```
ancify [-h] [-v] {init,project,call,evaluate,run} ...
```

## Global options

| Flag | Description |
|------|-------------|
| `-h`, `--help` | Show help and exit |
| `-v`, `--verbose` | Enable debug-level logging |

## Subcommands

### `ancify init`

Generate a template configuration file.

```bash
ancify init                    # writes config.yaml
ancify init -o my_config.yaml  # writes to custom path
```

| Option | Description |
|--------|-------------|
| `-o`, `--output` | Output path (default: `config.yaml`) |

### `ancify project`

Run **Phase 1**: project outgroup alignments onto focal-species coordinates.

```bash
ancify project -c config.yaml
ancify project -c config.yaml -n 8    # use 8 workers
```

Creates `<work_dir>/projected/<species>/<chrom>.fa` for every outgroup species
and chromosome.

| Option | Description |
|--------|-------------|
| `-c`, `--config` | Path to YAML config file (required) |
| `-n`, `--num-cpus` | Override `num_cpus` from config |

### `ancify call`

Run **Phase 2**: infer ancestral alleles from projected sequences.

```bash
ancify call -c config.yaml
```

Reads projected FASTA files and writes ancestral FASTA files to `<output_dir>/`.

| Option | Description |
|--------|-------------|
| `-c`, `--config` | Path to YAML config file (required) |
| `-n`, `--num-cpus` | Override `num_cpus` from config |

### `ancify evaluate`

Run **Phase 3**: evaluate ancestral calls against a reference and/or VCF data.

```bash
ancify evaluate -c config.yaml
```

Requires the `evaluation` block in the config and the `scikit-allel` package
for VCF comparison. Writes per-chromosome summary files to
`<output_dir>/evaluation/`.

| Option | Description |
|--------|-------------|
| `-c`, `--config` | Path to YAML config file (required) |
| `-n`, `--num-cpus` | Override `num_cpus` from config |

### `ancify run`

Run **all phases** end-to-end (project → call → evaluate).

```bash
ancify run -c config.yaml
ancify -v run -c config.yaml -n 4    # verbose, 4 workers
```

Phase 3 is only run if the `evaluation` block is present in the config.

| Option | Description |
|--------|-------------|
| `-c`, `--config` | Path to YAML config file (required) |
| `-n`, `--num-cpus` | Override `num_cpus` from config |

## Examples

Run the full pipeline:

```bash
ancify run -c example_configs/hg38_bcgm.yaml
```

Run only projection, then call separately:

```bash
ancify project -c config.yaml
ancify call -c config.yaml
```

Debug a failing run:

```bash
ancify -v run -c config.yaml -n 1 2>&1 | tee ancify.log
```
