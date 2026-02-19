# Configuration Reference

All pipeline behaviour is controlled by a single YAML configuration file.
Generate a starter template with:

```bash
ancify init -o config.yaml
```

## Complete annotated config

```yaml
# Informal label for the focal species (used in log messages).
focal_species: human

# Tab-separated file: chrom_name <TAB> length [<TAB> ...]
# Additional columns beyond the first two are ignored.
chromosome_lengths: chromoLens.txt

# Optional: restrict to specific chromosomes.
# If omitted, every entry in the lengths file is processed.
chromosomes:
  - chr1
  - chr2
  - chrX

outgroups:
  # Inner outgroup: closely related species.
  # A majority vote among these determines the inner consensus.
  inner:
    - name: bonobo
      alignment: hg38.panPan3.net.axt.gz
    - name: chimp
      alignment: hg38.panTro6.net.axt.gz
    - name: gorilla
      alignment: hg38.gorGor6.net.axt.gz

  # Outer outgroup: distantly related species.
  # An independent check on the inner consensus.
  outer:
    - name: macaque
      alignment: hg38.rheMac10.net.axt.gz

# Working directory for intermediate projected sequences.
# Projected files are stored in <work_dir>/projected/<species>/<chrom>.fa
work_dir: .

# Output directory for final ancestral FASTA files.
output_dir: ./ancestral_calls

# Minimum allele count to accept a majority-vote consensus.
# With 3 inner species: 1 = any single species suffices,
# 2 = at least 2 must agree, 3 = all must agree.
min_inner_freq: 1
min_outer_freq: 1

# Number of parallel worker processes.
num_cpus: 24

# Optional evaluation block (Phase 3).
evaluation:
  # Directory containing reference ancestral FASTA files.
  reference_dir: ./ensembl_ancestor/
  # Filename pattern. Supports {chrom} and {chrom_id} placeholders.
  reference_pattern: "homo_sapiens_ancestor_{chrom_id}.fa"

  # Directory containing VCF files.
  vcf_dir: ./vcf/
  vcf_pattern: "ALL.chr{chrom_id}.vcf.gz"
```

## Field reference

### Required fields

| Field | Type | Description |
|-------|------|-------------|
| `focal_species` | string | Label for the focal species (cosmetic only) |
| `chromosome_lengths` | path | Tab-separated file with at least two columns: chromosome name and length |
| `outgroups.inner` | list | One or more closely related outgroup species |
| `outgroups.outer` | list | One or more distantly related outgroup species |

Each outgroup entry requires:

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Species identifier (used in directory names) |
| `alignment` | path | Path to the net AXT pairwise alignment file |

### Optional fields

| Field | Default | Description |
|-------|---------|-------------|
| `chromosomes` | all from lengths file | List of chromosomes to process |
| `work_dir` | `.` | Working directory for projected sequences |
| `output_dir` | `./ancestral_calls` | Output directory for ancestral FASTAs |
| `min_inner_freq` | `1` | Minimum allele count for inner majority vote |
| `min_outer_freq` | `1` | Minimum allele count for outer majority vote |
| `num_cpus` | `4` | Number of parallel worker processes |
| `evaluation` | none | Optional evaluation settings (see below) |

### Evaluation fields

All evaluation fields are optional. If the `evaluation` block is omitted, Phase 3 is skipped.

| Field | Default | Description |
|-------|---------|-------------|
| `reference_dir` | none | Directory with reference ancestral FASTA files |
| `reference_pattern` | `{chrom}.fa` | Filename pattern for reference files |
| `vcf_dir` | none | Directory with VCF files |
| `vcf_pattern` | `{chrom}.vcf.gz` | Filename pattern for VCF files |

### Pattern placeholders

Evaluation filename patterns support two placeholders:

- **`{chrom}`** — the full chromosome name as it appears in the config (e.g., `chr1`, `2L`).
- **`{chrom_id}`** — the name with any leading `chr` prefix stripped (e.g., `1`, `X`).

This is useful when the ancestral reference uses a different naming convention than your focal genome. For example, ancify uses `chr1` but Ensembl files are named `homo_sapiens_ancestor_1.fa`, so you'd use:

```yaml
reference_pattern: "homo_sapiens_ancestor_{chrom_id}.fa"
```

## The chromosome lengths file

This is a simple tab-separated text file. Only the first two columns matter:

```
chr1    248956422
chr2    242193529
chrX    156040895
```

Additional columns (GenBank accession, RefSeq ID, etc.) are silently ignored.

### How to create one

From a FASTA index:

```bash
samtools faidx reference.fa
cut -f1,2 reference.fa.fai > chromoLens.txt
```

From UCSC:

```bash
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A \
  -e "SELECT chrom, size FROM chromInfo" hg38 > chromoLens.txt
```

## Validation

When a config is loaded, ancify validates:

- At least one inner outgroup species is specified.
- At least one outer outgroup species is specified.
- The chromosome lengths file exists.
- All alignment files exist.

If any check fails, a clear error message is printed before any processing begins.
