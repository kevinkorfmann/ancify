# Scripts

Test ancify end-to-end on real UCSC data. One unified script handles downloading,
config generation, and running; thin wrappers set method/backend for convenience.

## Layout

```
scripts/
├── run_hg38.sh             # unified runner (all env knobs)
├── compare_methods.sh      # run voting + parsimony + likelihood, compare output
├── visualize_comparison.py # generate figures from comparison output
├── README.md
└── examples/
    ├── chr22_voting.sh         # quick test: chr22, voting, CPU
    ├── chr22_voting_gpu.sh     # quick test: chr22, voting, GPU
    ├── chr22_parsimony.sh      # chr22, Fitch parsimony
    ├── chr22_likelihood.sh     # chr22, likelihood (HKY85)
    ├── chr22_ml.sh             # chr22, LightGBM (needs ML_MODEL_PATH)
    ├── full_voting.sh          # full genome, voting
    ├── full_voting_gpu.sh      # full genome, voting, GPU
    ├── full_parsimony.sh       # full genome, parsimony
    └── full_likelihood.sh      # full genome, likelihood
```

## Prerequisites

- `ancify` installed (`pip install .`)
- `wget` and `curl`
- For GPU: PyTorch with CUDA
- For ML: trained model (`ancify train -c config.yaml -o model.lgb`)

## Quick start

```bash
# Default: chr22, voting, auto backend
./scripts/run_hg38.sh

# Compare all three methods on chr22 (subdirs for easy diff)
./scripts/compare_methods.sh

# Or use a wrapper
./scripts/examples/chr22_parsimony.sh
```

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `CHROM` | `chr22` | Chromosome to process, or `all` for chr1-22 + chrX |
| `METHOD` | `voting` | `voting`, `parsimony`, `likelihood`, or `ml` |
| `BACKEND` | `auto` | `auto`, `cpu`, or `gpu` |
| `ANCIFY_CPUS` | `4` | Number of parallel workers |
| `WORK_DIR` | `ancify_test` | Working directory (downloads + output) |
| `ML_MODEL_PATH` | — | Path to `.lgb` model (required for `ml`) |

## Examples

```bash
# GPU + likelihood on chr22
BACKEND=gpu METHOD=likelihood ./scripts/run_hg38.sh

# Full genome on GPU with 24 cores
CHROM=all BACKEND=gpu ANCIFY_CPUS=24 ./scripts/run_hg38.sh

# ML on chr22
ML_MODEL_PATH=model.lgb ./scripts/examples/chr22_ml.sh

# Compare methods on a different chromosome
CHROM=chr1 ./scripts/compare_methods.sh
```

## Output structure

```
ancify_compare/output_chr22/
├── voting/chr22.fa
├── parsimony/chr22.fa
├── likelihood/chr22.fa
└── figures/
    ├── confidence_breakdown.png    # stacked bar: high/low/unresolved/missing
    ├── pairwise_agreement.png      # heatmap of agreement rates
    ├── disagreement_windows.png    # sliding-window disagreement along chrom
    └── called_site_overlap.png     # which methods call which sites (3 methods)
```

`compare_methods.sh` prints per-method coverage stats and pairwise agreement rates,
then generates the figures above (requires `matplotlib`: `pip install 'ancify[evaluate]'`).

You can also run the visualization separately on existing output:

```bash
python scripts/visualize_comparison.py ancify_compare/output_chr22 chr22 voting parsimony likelihood
```
