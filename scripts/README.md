# Scripts

Test ancify end-to-end on real data. Each species has its own runner script that
handles downloading (where available), config generation, and pipeline execution.
Examples are organised by species, each running **one chromosome** for a quick test.

## Layout

```
scripts/
├── run_hg38.sh              # human (hg38) — downloads from UCSC
├── run_mm39.sh              # mouse (mm39) — downloads from UCSC
├── run_dm6.sh               # Drosophila melanogaster (dm6) — downloads from UCSC
├── run_brassica.sh          # Brassica rapa — requires pre-downloaded data
├── compare_methods.sh       # run voting + parsimony + likelihood, compare output
├── visualize_comparison.py  # generate figures from comparison output
├── README.md
└── examples/
    ├── human/
    │   ├── run.sh               # chr22, voting, CPU
    │   ├── run_parsimony.sh     # chr22, Fitch parsimony
    │   ├── run_likelihood.sh    # chr22, likelihood (HKY85)
    │   └── run_voting_gpu.sh    # chr22, voting, GPU
    ├── mouse/
    │   └── run.sh               # chr19, voting, CPU
    ├── drosophila/
    │   └── run.sh               # chrom 4, voting, CPU
    └── brassica_rapa/
        └── run.sh               # A01, voting, CPU (needs data)
```

## Install (uv or pip)

From the repo root:

```bash
# With uv (recommended)
uv sync

# Or with pip
pip install .
```

For visualization (after `compare_methods.sh`): `uv sync --extra evaluate` or `pip install 'ancify[evaluate]'`.
For GPU: install PyTorch with CUDA (e.g. `pip install torch` with a CUDA-enabled build).
For ML: `uv sync --extra ml` or `pip install 'ancify[ml]'`.

## Prerequisites

- `ancify` installed (see above)
- `wget` and `curl`
- For GPU: PyTorch with CUDA
- For ML: trained model (`ancify train -c config.yaml -o model.lgb`)
- For Brassica rapa: alignment files placed in `WORK_DIR` (see `run_brassica.sh` header)

## Quick start

```bash
# Human — chr22, voting, CPU
./scripts/examples/human/run.sh

# Mouse — chr19, voting, CPU
./scripts/examples/mouse/run.sh

# Drosophila — chrom 4, voting, CPU
./scripts/examples/drosophila/run.sh

# Brassica rapa — A01, voting (needs data in WORK_DIR first)
./scripts/examples/brassica_rapa/run.sh
```

Or call a runner directly with env overrides:

```bash
# Human chr22, parsimony
CHROM=chr22 METHOD=parsimony ./scripts/run_hg38.sh

# Mouse chr1 instead of chr19
CHROM=chr1 ./scripts/run_mm39.sh

# Drosophila 2L instead of 4
CHROM=2L ./scripts/run_dm6.sh
```

## Compare methods (human chr22)

```bash
# Voting + parsimony + likelihood on chr22, then figures
./scripts/compare_methods.sh
```

## Voting CPU vs GPU benchmark (human)

Compare run-time for **voting only** on several chromosomes, CPU vs GPU, and plot results:

```bash
# Default: chr20, chr21, chr22 (voting only)
./scripts/benchmark_voting_cpu_gpu.sh
```

Custom chromosomes:

```bash
CHROMOSOMES="chr21 chr22" ./scripts/benchmark_voting_cpu_gpu.sh
```

Timings are stratified by phase: **Phase 1** (project: AXT → projected FASTA) and **Phase 2** (voting: ancestral calling). Results go to `ancify_voting_bench/voting_timings.csv` (columns: chromosome, backend, phase1_sec, phase2_sec, time_sec) and a two-panel bar chart to `ancify_voting_bench/voting_cpu_vs_gpu.png`. Set `SKIP_PLOT=1` to skip the figure.

## GPU examples

```bash
# Human chr22, GPU
./scripts/examples/human/run_voting_gpu.sh

# Or any species with env override
BACKEND=gpu ./scripts/examples/mouse/run.sh
```

## Species runners — environment variables

All four runners share the same env-var interface:

| Variable | Default | Description |
|----------|---------|-------------|
| `CHROM` | species-specific | Chromosome to process (single) |
| `METHOD` | `voting` | `voting`, `parsimony`, `likelihood`, or `ml` |
| `BACKEND` | `auto` | `auto`, `cpu`, or `gpu` |
| `ANCIFY_CPUS` | `4` | Number of parallel workers |
| `WORK_DIR` | species-specific | Working directory (downloads + output) |
| `ML_MODEL_PATH` | — | Path to `.lgb` model (required for `ml`) |

Default chromosomes per species:

| Runner | Species | Default `CHROM` | Default `WORK_DIR` |
|--------|---------|-----------------|---------------------|
| `run_hg38.sh` | Human | `chr22` | `ancify_test` |
| `run_mm39.sh` | Mouse | `chr19` | `ancify_test_mouse` |
| `run_dm6.sh` | Drosophila | `4` | `ancify_test_dm6` |
| `run_brassica.sh` | Brassica rapa | `A01` | `ancify_test_brassica` |

## Output structure

Each run produces output under `WORK_DIR/output_<CHROM>/<METHOD>/`:

```
ancify_test_mouse/output_chr19/voting/chr19.fa
ancify_test_dm6/output_4/voting/4.fa
```

For `compare_methods.sh` (human chr22):

```
ancify_compare/output_chr22/
├── voting/chr22.fa
├── parsimony/chr22.fa
├── likelihood/chr22.fa
└── figures/
    ├── confidence_breakdown.png
    ├── pairwise_agreement.png
    ├── disagreement_windows.png
    └── called_site_overlap.png
```

You can also run the visualization separately on existing output:

```bash
python scripts/visualize_comparison.py ancify_compare/output_chr22 chr22 voting parsimony likelihood
```
