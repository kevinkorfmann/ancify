# Configuration Reference

Everything in ancify is controlled by a single YAML file. This page is the complete reference.

## Generate a starter config

```bash
ancify init -o config.yaml
```

This creates a fully annotated template. You can also start from one of the example configs in `example_configs/`.

---

## Complete annotated config

```yaml
focal_species: human

chromosome_lengths: chromoLens.txt

# Optional: restrict to specific chromosomes.
# If omitted, every entry in the lengths file is processed.
chromosomes:
  - chr1
  - chr2
  - chrX

outgroups:
  inner:
    - name: bonobo
      alignment: hg38.panPan3.net.axt.gz
    - name: chimp
      alignment: hg38.panTro6.net.axt.gz
    - name: gorilla
      alignment: hg38.gorGor6.net.axt.gz
  outer:
    - name: macaque
      alignment: hg38.rheMac10.net.axt.gz

work_dir: .
output_dir: ./ancestral_calls

min_inner_freq: 1
min_outer_freq: 1

num_cpus: 24

# Inference method: "voting" (default), "parsimony", or "ml".
method: voting

# --- Parsimony options (only needed when method: parsimony) ---
# tree: "(((bonobo,chimp),gorilla),macaque)"   # inline Newick, or:
# tree: species_tree.nwk                        # path to .nwk file

# --- Likelihood options (only needed when method: likelihood) ---
# substitution_model: HKY85         # JC69 (default), K80, HKY85, or GTR
# model_kappa: 2.0                  # transition/transversion ratio (K80, HKY85)
# model_base_freqs: [0.3, 0.2, 0.2, 0.3]  # equilibrium freqs (HKY85, GTR)
# model_rates: [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # exchangeability rates (GTR)
# likelihood_high_threshold: 0.8    # min posterior → uppercase (high conf)
# likelihood_low_threshold: 0.5     # min posterior → lowercase (low conf)

# --- ML options (only needed when method: ml) ---
# ml_model_path: model.lgb          # path to trained LightGBM model
# ml_training_reference: ./ensembl_ancestor/  # optional: supervised labels
# ml_high_threshold: 0.8            # min probability → uppercase (high conf)
# ml_low_threshold: 0.5             # min probability → lowercase (low conf)

# Performance backend: "auto", "cpu", or "gpu".
# "auto" uses GPU when PyTorch + CUDA is available, otherwise CPU.
backend: auto

# Optional: restrict to specific GPUs (default: use all).
# Example: gpu_devices: [0, 1, 2]

# Optional evaluation block (Phase 3).
evaluation:
  reference_dir: ./ensembl_ancestor/
  reference_pattern: "homo_sapiens_ancestor_{chrom_id}.fa"
  vcf_dir: ./vcf/
  vcf_pattern: "ALL.chr{chrom_id}.vcf.gz"
```

---

## Field reference

### Required fields

| Field | Type | Description |
|-------|------|-------------|
| `focal_species` | string | Label for the focal species (used in log messages only) |
| `chromosome_lengths` | path | Tab-separated file with at least two columns: chromosome name and length |
| `outgroups.inner` | list | One or more closely related outgroup species |
| `outgroups.outer` | list | One or more distantly related outgroup species |

Each outgroup entry requires:

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Species identifier (used in directory names for projected files) |
| `alignment` | path | Path to the net AXT pairwise alignment file (may be gzipped) |

### Optional fields

| Field | Default | Description |
|-------|---------|-------------|
| `chromosomes` | all from lengths file | List of chromosomes to process |
| `work_dir` | `.` | Working directory for intermediate projected sequences |
| `output_dir` | `./ancestral_calls` | Output directory for final ancestral FASTAs |
| `min_inner_freq` | `1` | Minimum allele count for inner majority vote |
| `min_outer_freq` | `1` | Minimum allele count for outer majority vote |
| `num_cpus` | `4` | Number of parallel worker processes |
| `backend` | `auto` | Compute backend: `"auto"`, `"cpu"`, or `"gpu"` (see {doc}`performance`) |
| `gpu_devices` | all available | List of GPU device IDs to use, e.g. `[0, 1, 2]` |
| `method` | `voting` | Ancestral inference method: `"voting"`, `"parsimony"`, `"likelihood"`, or `"ml"` (see {doc}`algorithm`) |
| `tree` | none | Newick tree of outgroup species (required when `method: parsimony` or `likelihood`). Inline string or path to `.nwk` file. |
| `substitution_model` | `JC69` | Substitution model for likelihood method: `"JC69"`, `"K80"`, `"HKY85"`, or `"GTR"`. |
| `model_kappa` | `2.0` | Transition/transversion ratio κ (used by K80 and HKY85). |
| `model_base_freqs` | uniform | Equilibrium base frequencies `[π_A, π_C, π_G, π_T]` (used by HKY85 and GTR). Must sum to 1. |
| `model_rates` | none | Six GTR exchangeability rates `[AC, AG, AT, CG, CT, GT]` (required when `substitution_model: GTR`). |
| `likelihood_high_threshold` | `0.8` | Minimum posterior probability for high (uppercase) confidence (likelihood method). |
| `likelihood_low_threshold` | `0.5` | Minimum posterior probability for low (lowercase) confidence (likelihood method). Below this → `n`. |
| `ml_model_path` | none | Path to a trained LightGBM model file (required when `method: ml`). Produced by `ancify train`. |
| `ml_training_reference` | none | Directory of reference ancestral FASTAs for supervised training (optional; used by `ancify train`). |
| `ml_high_threshold` | `0.8` | Minimum predicted probability for high (uppercase) confidence. |
| `ml_low_threshold` | `0.5` | Minimum predicted probability for low (lowercase) confidence. Below this → `n`. |
| `evaluation` | none | Optional evaluation block (see below) |

### Evaluation fields

All evaluation fields are optional. If the `evaluation` block is omitted entirely, Phase 3 is skipped.

| Field | Default | Description |
|-------|---------|-------------|
| `evaluation.reference_dir` | none | Directory with reference ancestral FASTA files |
| `evaluation.reference_pattern` | `{chrom}.fa` | Filename pattern for reference files |
| `evaluation.vcf_dir` | none | Directory with VCF files |
| `evaluation.vcf_pattern` | `{chrom}.vcf.gz` | Filename pattern for VCF files |

---

## Understanding key parameters

### `min_inner_freq`: the stringency dial

This is the single most important tuning parameter. It controls how many inner outgroup species must agree before a consensus is accepted.

```text
  With 3 inner species (e.g. bonobo, chimp, gorilla):

  min_inner_freq=1     min_inner_freq=2     min_inner_freq=3
  ──────────────────   ──────────────────   ──────────────────
  Any 1 species        At least 2 of 3      All 3 must agree
  suffices             must agree

  ✓ Maximum coverage   ✓ Balanced            ✓ Maximum stringency
  ✗ Lower accuracy     ✓ Good for most       ✗ Lower coverage
                         use cases            ✓ Highest accuracy
```

**Recommendation:** Start with the default (`1`) and increase if you see higher-than-expected disagreement rates in the evaluation.

### `num_cpus`: parallelism

Each chromosome is processed independently. `num_cpus` controls how many chromosomes are processed simultaneously.

```text
  Memory usage ≈ num_cpus × (num_outgroups) × (avg_chrom_length)

  Example: 24 CPUs × 4 outgroups × 150 Mb avg = ~14.4 GB peak
```

If you are running on a machine with limited memory, reduce `num_cpus`. For Phase 1 (projection), memory is modest; the bottleneck is Phase 2 (calling), which loads all projected sequences for each active chromosome.

### `method`: choosing your inference approach

ancify supports four methods for inferring the ancestral allele at each position.
You switch between them with a single line in your YAML:

```yaml
method: voting      # default
```

```text
  ┌──────────────┬────────────────────────────────────────────────────────────┐
  │ method       │ Description                                                │
  ├──────────────┼────────────────────────────────────────────────────────────┤
  │ voting       │ Two-tier majority vote. Inner outgroups vote first; outer  │
  │ (default)    │ provides independent check. Fast, no extra deps.           │
  ├──────────────┼────────────────────────────────────────────────────────────┤
  │ parsimony    │ Fitch (1971) algorithm on a Newick phylogenetic tree.      │
  │              │ Resolves many "unresolved" (n) sites that voting misses.   │
  │              │ Requires: tree field.                                      │
  ├──────────────┼────────────────────────────────────────────────────────────┤
  │ likelihood   │ Felsenstein (1981) pruning with an explicit substitution   │
  │              │ model. Computes posterior probabilities at the root using   │
  │              │ branch lengths. Requires: tree with branch lengths, scipy. │
  ├──────────────┼────────────────────────────────────────────────────────────┤
  │ ml           │ LightGBM gradient-boosted classifier. Learns substitution  │
  │              │ biases (CpG, GC context, etc.) from data. Yields calibrated│
  │              │ probability scores rather than binary high/low.            │
  │              │ Requires: pip install 'ancify[ml]', a trained model.       │
  └──────────────┴────────────────────────────────────────────────────────────┘
```

**Not sure which to use?** Start with `voting`. It requires no extra configuration
and gives >99.9% accuracy for most species. Switch to `parsimony` if you have a
reliable tree and want to resolve more "unresolved" sites. Use `likelihood` if you
have a tree with branch lengths and want calibrated posterior probabilities from
an explicit substitution model. Use `ml` if you have an external reference sequence
for training labels and want to learn context-dependent patterns.
See {doc}`algorithm` for a detailed comparison with worked examples.

#### Parsimony YAML

```yaml
method: parsimony
tree: "(((bonobo,chimp),gorilla),macaque)"
```

Leaf names in the Newick string must exactly match the `name` fields of your
outgroup entries. You can also point to a file:

```yaml
tree: species_tree.nwk
```

#### Likelihood YAML

```yaml
method: likelihood
tree: "(((bonobo:0.008,chimp:0.008):0.002,gorilla:0.009):0.020,macaque:0.038)"

# Substitution model (default: JC69)
substitution_model: HKY85
model_kappa: 2.0
model_base_freqs: [0.3, 0.2, 0.2, 0.3]

# Optional: tune confidence thresholds (defaults shown)
likelihood_high_threshold: 0.8   # posterior ≥ 0.8 → UPPERCASE (high conf)
likelihood_low_threshold: 0.5    # posterior ≥ 0.5 → lowercase (low conf)
                                 # below 0.5 → n (unresolved)
```

Branch lengths are required — they control how much influence each leaf has.
Trees from RAxML, IQ-TREE, or similar tools provide lengths in substitutions
per site, which is exactly what the models expect. For GTR, also supply six
exchangeability rates:

```yaml
substitution_model: GTR
model_rates: [1.0, 2.0, 1.0, 1.0, 2.0, 1.0]   # [AC, AG, AT, CG, CT, GT]
model_base_freqs: [0.3, 0.2, 0.2, 0.3]
```

#### ML YAML (two-step workflow)

**Step 1: Train the model** (run once, before calling):

```yaml
# Add to config.yaml — optionally point to a reference for supervised labels:
ml_training_reference: /data/ensembl_ancestor/   # omit for self-supervised
```

```bash
ancify train -c config.yaml -o model.lgb
```

**Step 2: Call with the trained model**:

```yaml
method: ml
ml_model_path: model.lgb

# Optional: tune confidence thresholds (defaults shown)
ml_high_threshold: 0.8   # predicted probability ≥ 0.8 → UPPERCASE (high conf)
ml_low_threshold: 0.5    # predicted probability ≥ 0.5 → lowercase (low conf)
                         # below 0.5 → n (unresolved)
```

```bash
ancify call -c config.yaml
```

Install ML dependencies first if you haven't:

```bash
pip install 'ancify[ml]'
```

### `backend`: the compute engine

This field selects which execution backend ancify uses for the heavy computation
in Phases 1 and 2. See {doc}`performance` for the full details.

```text
  backend: "auto"             backend: "cpu"             backend: "gpu"
  ──────────────────          ──────────────             ──────────────
  Uses GPU if PyTorch         Forces CPU-only            Forces GPU mode
  + CUDA is detected,         (vectorized NumPy).        (requires PyTorch
  else falls back to          Still much faster than     + CUDA).
  vectorized CPU.             the unvectorized path.
```

When `backend: "gpu"` (or `"auto"` with CUDA detected) is active, chromosomes
are distributed round-robin across all available GPUs. Use `gpu_devices` to
restrict which ones are used:

```yaml
# Use all available GPUs (default)
backend: auto

# Use specific GPUs only
backend: auto
gpu_devices: [0, 1, 2]   # three A100s

# Single GPU, explicit
backend: gpu
gpu_devices: [0]

# Force CPU-only even if a GPU is present
backend: cpu
```

Install PyTorch with CUDA to unlock the GPU path:

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu128
```

See {doc}`performance` for supported hardware, memory budgets, and benchmarks.

### Choosing chromosomes

By default, ancify processes every chromosome in the lengths file. This includes autosomes, sex chromosomes, and scaffolds. To process only specific chromosomes:

```yaml
chromosomes:
  - chr1
  - chr2
  - chrX
```

:::{tip}
For an initial test run, try processing just one small chromosome (e.g. `chr22` for humans or `chr19` for mouse) to verify everything works before committing to a full-genome run.
:::

---

## Pattern placeholders

Evaluation filename patterns support two placeholders:

| Placeholder | Example input: `chr1` | Example input: `2L` |
|-------------|----------------------|---------------------|
| `{chrom}` | `chr1` | `2L` |
| `{chrom_id}` | `1` (strips `chr` prefix) | `2L` (no `chr` prefix to strip) |

This handles the common case where your focal genome uses `chr1` but reference files use `1`:

```yaml
reference_pattern: "homo_sapiens_ancestor_{chrom_id}.fa"
#                                          ↑ becomes "1" for chr1
```

---

## The chromosome lengths file

A simple tab-separated text file. Only the first two columns matter:

```text
chr1    248956422
chr2    242193529
chrX    156040895
```

Additional columns (GenBank accession, RefSeq ID, etc.) are silently ignored.

### How to create one

**From a FASTA index:**

```bash
samtools faidx reference.fa
cut -f1,2 reference.fa.fai > chromoLens.txt
```

**From UCSC MySQL:**

```bash
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A \
  -e "SELECT chrom, size FROM chromInfo" hg38 > chromoLens.txt
```

**From an NCBI assembly report:**

```bash
grep -v "^#" assembly_report.txt | awk -F'\t' '{print $1, $9}' > chromoLens.txt
```

---

## Validation

When a config is loaded, ancify validates:

- At least one inner outgroup species is specified (voting, ml).
- At least one outer outgroup species is specified (voting, ml).
- At least one outgroup species is specified (parsimony, likelihood).
- The chromosome lengths file exists and is readable.
- All alignment files exist and are accessible.
- If `method` is `parsimony` or `likelihood`: a `tree` field is present and leaf names match outgroup names.
- If `method` is `likelihood`: the substitution model name is valid, GTR has 6 rates, base frequencies (if provided) have 4 entries summing to ~1, and thresholds satisfy 0 ≤ low ≤ high ≤ 1.

If any check fails, a clear error message is printed before any processing begins. This prevents wasting hours on Phase 1 only to discover a typo in a file path.

---

## Config recipes

### Minimal (2 species)

The simplest possible configuration — one inner, one outer:

```yaml
focal_species: mouse
chromosome_lengths: mm39.chromLens.txt
outgroups:
  inner:
    - name: rat
      alignment: mm39.rn7.net.axt.gz
  outer:
    - name: rabbit
      alignment: mm39.oryCun2.net.axt.gz
output_dir: ./mouse_ancestral
```

### Maximal (many species, strict settings, evaluation)

```yaml
focal_species: human
chromosome_lengths: chromoLens.txt
chromosomes: [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9,
              chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17,
              chr18, chr19, chr20, chr21, chr22, chrX]
outgroups:
  inner:
    - name: bonobo
      alignment: hg38.panPan3.net.axt.gz
    - name: chimp
      alignment: hg38.panTro6.net.axt.gz
    - name: gorilla
      alignment: hg38.gorGor6.net.axt.gz
  outer:
    - name: macaque
      alignment: hg38.rheMac10.net.axt.gz
work_dir: /scratch/polarization
output_dir: /results/human_ancestral
min_inner_freq: 2
min_outer_freq: 1
num_cpus: 32
backend: auto
gpu_devices: [0, 1, 2]
evaluation:
  reference_dir: /data/ensembl_ancestor/
  reference_pattern: "homo_sapiens_ancestor_{chrom_id}.fa"
  vcf_dir: /data/1kg/
  vcf_pattern: "ALL.chr{chrom_id}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
```

### Fitch parsimony

The same outgroup configuration with the phylogenetic tree-based method:

```yaml
focal_species: human
chromosome_lengths: chromoLens.txt
outgroups:
  inner:
    - name: bonobo
      alignment: hg38.panPan3.net.axt.gz
    - name: chimp
      alignment: hg38.panTro6.net.axt.gz
    - name: gorilla
      alignment: hg38.gorGor6.net.axt.gz
  outer:
    - name: macaque
      alignment: hg38.rheMac10.net.axt.gz
method: parsimony
tree: "(((bonobo,chimp),gorilla),macaque)"
output_dir: ./human_ancestral_parsimony
```

Leaf names in the Newick tree must match the `name` fields of your outgroup entries. When `method` is `parsimony`, the inner/outer distinction is ignored — all outgroups are placed on the tree and the Fitch algorithm determines their relative weights.

You can also point `tree` to an external file:

```yaml
tree: species_tree.nwk
```

### Likelihood

The same outgroup configuration with the likelihood method and an HKY85 substitution model:

```yaml
focal_species: human
chromosome_lengths: chromoLens.txt
outgroups:
  inner:
    - name: bonobo
      alignment: hg38.panPan3.net.axt.gz
    - name: chimp
      alignment: hg38.panTro6.net.axt.gz
    - name: gorilla
      alignment: hg38.gorGor6.net.axt.gz
  outer:
    - name: macaque
      alignment: hg38.rheMac10.net.axt.gz
method: likelihood
tree: "(((bonobo:0.008,chimp:0.008):0.002,gorilla:0.009):0.020,macaque:0.038)"
substitution_model: HKY85
model_kappa: 2.0
model_base_freqs: [0.3, 0.2, 0.2, 0.3]
output_dir: ./human_ancestral_likelihood
```

When `method` is `likelihood`, branch lengths in the tree are critical — they determine
how much weight each leaf's observation receives. The inner/outer distinction is ignored;
all outgroups are placed on the tree.

### ML classifier

The ML method requires a two-step workflow: train first, then call.

**Step 1: Train** (uses high-confidence voting sites as labels by default):

```yaml
focal_species: human
chromosome_lengths: chromoLens.txt
outgroups:
  inner:
    - name: bonobo
      alignment: hg38.panPan3.net.axt.gz
    - name: chimp
      alignment: hg38.panTro6.net.axt.gz
    - name: gorilla
      alignment: hg38.gorGor6.net.axt.gz
  outer:
    - name: macaque
      alignment: hg38.rheMac10.net.axt.gz
output_dir: ./human_ancestral_ml
```

```bash
ancify train -c config.yaml -o model.lgb
```

To use an external reference (e.g. Ensembl EPO) as training labels instead, add:

```yaml
ml_training_reference: /data/ensembl_ancestor/
```

The reference directory should contain one FASTA per chromosome (e.g. `chr1.fa`).

**Step 2: Call** (add the trained model path and set `method: ml`):

```yaml
method: ml
ml_model_path: model.lgb

# Optional: tune confidence thresholds (defaults shown)
ml_high_threshold: 0.8
ml_low_threshold: 0.5
```

```bash
ancify call -c config.yaml
```

### Quick test run (single chromosome)

```yaml
focal_species: human
chromosome_lengths: chromoLens.txt
chromosomes:
  - chr22
outgroups:
  inner:
    - name: chimp
      alignment: hg38.panTro6.net.axt.gz
  outer:
    - name: macaque
      alignment: hg38.rheMac10.net.axt.gz
output_dir: ./test_run
num_cpus: 1
```
