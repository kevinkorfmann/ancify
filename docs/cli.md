# Command-Line Interface

ancify provides a single CLI command with several subcommands. It can be invoked either as `ancify` (after installation) or as `python -m ancify`.

---

## Synopsis

```text
ancify [-h] [-v] {init,project,call,evaluate,run} ...
```

## Global options

| Flag | Description |
|------|-------------|
| `-h`, `--help` | Show help and exit |
| `-v`, `--verbose` | Enable debug-level logging (prints per-chromosome progress, timing, etc.) |

---

## Subcommands

### `ancify init` — Generate a config template

```bash
ancify init                    # writes config.yaml to current directory
ancify init -o my_config.yaml  # writes to a custom path
```

| Option | Description |
|--------|-------------|
| `-o`, `--output` | Output path (default: `config.yaml`) |

This creates a fully annotated YAML template with all fields documented. It is a good starting point for any new species.

### `ancify project` — Phase 1: Coordinate projection

```bash
ancify project -c config.yaml
ancify project -c config.yaml -n 8    # override num_cpus
```

| Option | Description |
|--------|-------------|
| `-c`, `--config` | Path to YAML config file (**required**) |
| `-n`, `--num-cpus` | Override `num_cpus` from config |

**What it does:** For each outgroup species and each chromosome, reads the pairwise alignment file and projects the outgroup's bases onto the focal genome's coordinate system.

**Output:** `<work_dir>/projected/<species>/<chrom>.fa` for every (species, chromosome) pair.

**Runtime:** This is the slow step (hours for a full human genome). Each AXT file is streamed sequentially per chromosome.

### `ancify call` — Phase 2: Ancestral state inference

```bash
ancify call -c config.yaml
ancify call -c config.yaml -n 4
```

| Option | Description |
|--------|-------------|
| `-c`, `--config` | Path to YAML config file (**required**) |
| `-n`, `--num-cpus` | Override `num_cpus` from config |

**What it does:** Reads projected FASTA files, computes two-tier majority voting at every position, and writes confidence-encoded ancestral FASTA files.

**Output:** `<output_dir>/<chrom>.fa` for each chromosome.

**Runtime:** Fast (minutes for a full human genome).

### `ancify evaluate` — Phase 3: Evaluation

```bash
ancify evaluate -c config.yaml
```

| Option | Description |
|--------|-------------|
| `-c`, `--config` | Path to YAML config file (**required**) |
| `-n`, `--num-cpus` | Override `num_cpus` from config |

**What it does:** Compares ancestral calls against a reference ancestral sequence and/or VCF variant data.

**Requires:** The `evaluation` block in the config, and `scikit-allel` for VCF comparison.

**Output:** Per-chromosome evaluation files in `<output_dir>/evaluation/`.

### `ancify run` — All phases end-to-end

```bash
ancify run -c config.yaml
ancify run -c config.yaml -n 24
ancify -v run -c config.yaml          # verbose output
```

| Option | Description |
|--------|-------------|
| `-c`, `--config` | Path to YAML config file (**required**) |
| `-n`, `--num-cpus` | Override `num_cpus` from config |

**What it does:** Runs Phase 1 → Phase 2 → Phase 3 in sequence. Phase 3 is skipped if the `evaluation` block is absent.

---

## Workflow patterns

### Standard full run

The most common usage — run everything from start to finish:

```bash
ancify run -c config.yaml
```

### Iterate on Phase 2 settings

Phase 1 is expensive. Once the projected files exist, you can re-run Phase 2 with different settings without redoing Phase 1:

```bash
# First time: run everything
ancify run -c config.yaml

# Later: tweak min_inner_freq and re-call
# (edit config.yaml to change min_inner_freq)
ancify call -c config.yaml
```

### Debug a failing run

```bash
ancify -v run -c config.yaml -n 1 2>&1 | tee ancify.log
```

- `-v` enables verbose logging (per-chromosome progress)
- `-n 1` uses a single worker (easier to read output, avoids interleaved logs)
- `tee` saves output to a file while still printing to screen

### Test with a single chromosome first

Edit your config to process only one small chromosome:

```yaml
chromosomes:
  - chr22    # smallest human autosome
```

Then run normally. This lets you verify the setup in minutes instead of hours.

### Process species independently (Phase 1)

For large genomes, you can parallelize Phase 1 across machines by running separate configs for each outgroup, then combining:

```bash
# On machine 1:
ancify project -c config_bonobo_only.yaml

# On machine 2:
ancify project -c config_chimp_only.yaml

# Then combine on one machine:
ancify call -c config_all.yaml
```

Phase 2 reads from whatever projected files exist in `<work_dir>/projected/`.

---

## Examples

Run the included human example:

```bash
ancify run -c example_configs/hg38_bcgm.yaml
```

Run only projection, then call separately:

```bash
ancify project -c config.yaml
ancify call -c config.yaml
```

Generate a config and inspect it:

```bash
ancify init -o my_species.yaml
cat my_species.yaml
```
