# ancify

**Infer ancestral alleles for any species using outgroup alignments.**

[![Documentation](https://img.shields.io/badge/docs-Read%20the%20Docs-blue)](https://ancify.readthedocs.io)
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen)](https://github.com/kevinkorfmann/ancify/actions)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)

ancify is a config-driven Python pipeline that determines the ancestral state at every position in a reference genome by comparing pairwise alignments from multiple outgroup species. It supports **four inference methods**: **two-tier inner/outer outgroup voting**, **Fitch parsimony** on a phylogenetic tree, **likelihood-based reconstruction** (Felsenstein pruning with substitution models JC69/K80/HKY85/GTR), and a **LightGBM classifier** — all with case-encoded confidence levels.

**Full documentation:** [ancify.readthedocs.io](https://ancify.readthedocs.io) — includes a population genetics background primer, step-by-step tutorials, algorithm deep dives, and a species adaptation guide.

---

## For the Impatient

```bash
pip install .
ancify init -o config.yaml   # generate template config
# edit config.yaml with your species, alignments, and paths
ancify run -c config.yaml    # run everything
```

That's it. Your ancestral FASTA files appear in the configured `output_dir`. Uppercase = high confidence, lowercase = low confidence.

### Quick Example (Human, hg38)

```bash
ancify run -c example_configs/hg38_bcgm.yaml
```

This polarizes the human genome using bonobo + chimp + gorilla (inner) and macaque (outer), producing one ancestral FASTA per chromosome.

---

## What It Does

```
Net AXT alignments          Projected sequences         Ancestral FASTA
(outgroup vs focal)   -->   (in focal coordinates)  --> (with confidence)

  hg38.panTro6.axt.gz       projected/chimp/chr1.fa
  hg38.panPan3.axt.gz  -->  projected/bonobo/chr1.fa --> ancestral/chr1.fa
  hg38.gorGor6.axt.gz       projected/gorilla/chr1.fa
  hg38.rheMac10.axt.gz      projected/macaque/chr1.fa

  Phase 1: project           Phase 2: call              Done.
```

**Phase 1** projects each outgroup alignment onto the focal genome's coordinates.
**Phase 2** infers the ancestral allele at every position using your chosen method (voting, parsimony, likelihood, or ML).
**Phase 3** (optional) evaluates calls against a reference and/or VCF variants.

### Confidence Encoding

**Voting method** (default):

| Character | Confidence | Meaning |
|-----------|-----------|---------|
| `ACGT` | High | Inner and outer outgroups agree |
| `acgt` | Low | Only one tier has data |
| `n` | Unresolved | Inner and outer disagree |
| `N` | Missing | No data from either tier |

**Parsimony method**:

| Character | Confidence | Meaning |
|-----------|-----------|---------|
| `ACGT` | High | Unique most-parsimonious root state |
| `acgt` | Low | Ambiguous root (multiple equally parsimonious states) |
| `N` | Missing | All outgroup leaves lack data |

**Likelihood method** (Felsenstein pruning):

| Character | Confidence | Meaning |
|-----------|-----------|---------|
| `ACGT` | High | Root posterior ≥ high threshold (default 0.8) |
| `acgt` | Low | Root posterior ≥ low threshold (default 0.5) |
| `n` | Unresolved | Root posterior &lt; low threshold |
| `N` | Missing | All outgroup leaves lack data |

**ML method**: same as voting (probability thresholds map to uppercase/lowercase/`n`/`N`).

---

## Installation

### Prerequisites

- **Python 3.8 or newer**
- **pip** (or another PEP 517–compatible installer)

Optional: a [virtual environment](https://docs.python.org/3/library/venv.html) is recommended to isolate dependencies.

### Core install (voting, parsimony, likelihood)

Clone the repository and install in editable or regular mode:

```bash
git clone https://github.com/kevinkorfmann/ancify.git
cd ancify
pip install .
```

This pulls in the core dependencies: **PyYAML**, **NumPy**, **SciPy** (for the likelihood method), and **LightGBM**. The **voting**, **parsimony**, and **likelihood** inference methods work out of the box.

**From a tarball or wheel (if you have a release):**

```bash
pip install ancify
# or
pip install dist/ancify-*.whl
```

### Optional extras

Install optional dependency groups with the `[extra]` syntax:

| Extra        | Purpose |
|-------------|---------|
| `evaluate` | Phase 3 evaluation: scikit-allel, matplotlib (reference & VCF comparison) |
| `fast`     | Faster gzip I/O via `isal` |
| `ml`       | ML method: LightGBM, scikit-learn (already in core; this adds sklearn explicitly) |
| `docs`     | Build Sphinx docs: sphinx, theme, myst-parser |
| `dev`      | Testing: pytest |
| `all`      | All of the above |

**Examples:**

```bash
# Core + evaluation (reference comparison, VCF concordance)
pip install '.[evaluate]'

# Core + faster gzip for large AXT files
pip install '.[fast]'

# Core + evaluation + fast I/O
pip install '.[evaluate,fast]'

# Everything (development and docs)
pip install '.[all]'
```

### GPU acceleration (voting method only)

The **voting** method can use PyTorch with CUDA to process full genomes in minutes instead of hours. The **parsimony**, **likelihood**, and **ML** methods run on CPU (multi-process).

1. Install ancify first (as above).
2. Install PyTorch with CUDA for your driver (example: CUDA 12.8):

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu128
```

3. In your config, use `backend: auto` (default) so ancify will use the GPU if available, or `backend: gpu` to require it.

See the [performance docs](https://ancify.readthedocs.io/en/latest/performance.html) for supported setups and benchmarks.

### Verify installation

```bash
ancify --help
python -c "from ancify.config import load_config; from ancify.likelihood import build_model; print('ok')"
```

---

### Installation quick reference

| What you want              | Command |
|----------------------------|---------|
| Minimal (voting, parsimony, likelihood) | `pip install .` |
| + Phase 3 evaluation       | `pip install '.[evaluate]'` |
| + Faster gzip              | `pip install '.[fast]'` |
| + GPU (voting only)        | `pip install .` then `pip install torch --index-url ...` |
| Dev/test                   | `pip install '.[dev]'` |

---

## Configuration

Everything is controlled by a single YAML file. Generate a starter template:

```bash
ancify init -o config.yaml
```

### Minimal Config

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

output_dir: ./ancestral_calls
num_cpus: 24
```

### Key Fields

| Field | Description |
|-------|-------------|
| `focal_species` | Label for the focal species (cosmetic) |
| `chromosome_lengths` | Tab-separated file: `chrom_name\tlength` |
| `chromosomes` | Optional list; defaults to all chroms in lengths file |
| `outgroups.inner` | Closely related species (majority vote) |
| `outgroups.outer` | Distantly related species (independent check) |
| `min_inner_freq` | Min count for inner majority vote (default: 1) |
| `min_outer_freq` | Min count for outer majority vote (default: 1) |
| `method` | `"voting"` (default), `"parsimony"`, `"likelihood"`, or `"ml"` |
| `tree` | Newick tree (required for parsimony/likelihood); use branch lengths for likelihood |
| `substitution_model` | For likelihood: `JC69`, `K80`, `HKY85`, or `GTR` (default: JC69) |
| `model_kappa`, `model_base_freqs`, `model_rates` | Likelihood model parameters (see docs) |
| `likelihood_high_threshold`, `likelihood_low_threshold` | Posterior thresholds for likelihood (defaults 0.8, 0.5) |
| `ml_model_path`, `ml_high_threshold`, `ml_low_threshold` | For ML method (train with `ancify train`) |
| `num_cpus` | Parallel workers (default: 4) |
| `evaluation` | Optional block for Phase 3 (reference + VCF comparison) |

---

## CLI Reference

```bash
ancify init     [-o FILE]           # generate template config
ancify project  -c CONFIG [-n N]    # Phase 1: project alignments
ancify call     -c CONFIG [-n N]    # Phase 2: call ancestral states
ancify evaluate -c CONFIG [-n N]    # Phase 3: evaluate (optional)
ancify run      -c CONFIG [-n N]    # all phases end-to-end
ancify train    -c CONFIG [-o MODEL] [-n N]   # train ML model (for method: ml)
```

Also: `python -m ancify run -c config.yaml`

---

## Works With Any Species

ancify is not tied to humans. It works with any focal species for which you have pairwise net AXT alignments (widely available from the [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html)).

### Mouse

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
```

### Drosophila

```yaml
focal_species: drosophila_melanogaster
chromosome_lengths: dm6.chromLens.txt
chromosomes: [2L, 2R, 3L, 3R, 4, X]
outgroups:
  inner:
    - name: simulans
      alignment: dm6.droSim2.net.axt.gz
    - name: sechellia
      alignment: dm6.droSec1.net.axt.gz
  outer:
    - name: yakuba
      alignment: dm6.droYak3.net.axt.gz
```

### Brassica rapa (plant)

```yaml
focal_species: brassica_rapa
chromosome_lengths: braRap1.chromLens.txt
outgroups:
  inner:
    - name: brassica_oleracea
      alignment: braRap1.braOleracea.net.axt.gz
  outer:
    - name: arabidopsis_thaliana
      alignment: braRap1.araTha1.net.axt.gz
```

See `example_configs/` for complete examples.

---

## How It Works

### Method 1: Two-tier voting (default)

1. **Inner consensus**: majority vote among closely related outgroup species (e.g. bonobo, chimp, gorilla).
2. **Outer consensus**: majority vote among distantly related outgroup species (e.g. macaque).
3. **Compare**:
   - Agree &rarr; **high confidence** (uppercase)
   - One missing &rarr; **low confidence** (lowercase, use the available call)
   - Disagree &rarr; **unresolved** (`n`)
   - Both missing &rarr; **missing** (`N`)

This two-tier approach guards against incomplete lineage sorting and lineage-specific substitutions. The outer outgroup provides an independent evolutionary check on the inner consensus.

### Method 2: Fitch parsimony

Instead of splitting outgroups into two tiers, you provide a **Newick phylogenetic tree** and ancify uses the **Fitch (1971) algorithm** to reconstruct the most parsimonious ancestral state at the root:

```yaml
method: parsimony
tree: "(((bonobo,chimp),gorilla),macaque)"
```

The tree topology determines how species are weighted, resolving cases that the voting method marks as unresolved. See the [algorithm docs](https://ancify.readthedocs.io/en/latest/algorithm.html) for a detailed walkthrough.

### Method 3: Likelihood (Felsenstein pruning)

With a **tree that has branch lengths**, you can use substitution-model-based reconstruction. ancify implements **Felsenstein’s pruning algorithm** and four models: **JC69**, **K80**, **HKY85**, and **GTR**. The root posterior probabilities are computed from the leaf data and branch lengths; confidence is encoded by thresholding the maximum posterior (e.g. ≥0.8 → uppercase, ≥0.5 → lowercase, &lt;0.5 → `n`).

```yaml
method: likelihood
tree: "(((bonobo:0.008,chimp:0.008):0.002,gorilla:0.009):0.020,macaque:0.038)"
substitution_model: HKY85
model_kappa: 2.0
```

Requires **SciPy** (included in the default install). See the [algorithm](https://ancify.readthedocs.io/en/latest/algorithm.html) and [configuration](https://ancify.readthedocs.io/en/latest/configuration.html) docs.

### Method 4: ML classifier

Train a LightGBM classifier on your data (or a reference), then use it for calling: `ancify train -c config.yaml -o model.lgb`, then set `method: ml` and `ml_model_path: model.lgb`. See the docs for feature description and calibration.

### Input: Net AXT Alignments

The pipeline reads pairwise **net AXT** alignment files from UCSC. These represent best-in-genome one-to-one alignments between the focal species and each outgroup. Download them from:

```
https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/vs<Outgroup>/
```

### Getting Chromosome Lengths

From UCSC:
```bash
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A \
  -e "SELECT chrom, size FROM chromInfo" hg38 > chromoLens.txt
```

Or from a FASTA index:
```bash
samtools faidx reference.fa
cut -f1,2 reference.fa.fai > chromoLens.txt
```

---

## Using the Output

### Look Up an Ancestral Allele

```python
from ancify.utils import read_fasta

_, seq = read_fasta("ancestral_calls/chr1.fa")
allele = seq[999999]  # 0-based index for position 1,000,000
print(f"Ancestral: {allele}, High confidence: {allele in 'ACGT'}")
```

### Polarize a VCF

```python
for variant in vcf:
    anc = seq[variant.POS - 1].upper()
    if anc == variant.REF:
        # REF is ancestral, ALT is derived
        ...
    elif anc == variant.ALT:
        # ALT is ancestral, REF is derived (flip frequencies)
        ...
```

### Python API

```python
from ancify.config import load_config
from ancify.project import run_projection
from ancify.ancestral import run_ancestral_calling, call_ancestral_base

# Run the full pipeline programmatically
cfg = load_config("config.yaml")
run_projection(cfg)
run_ancestral_calling(cfg)

# Or call the core function directly (voting)
base = call_ancestral_base(
    inner_bases=["A", "A", "G"],
    outer_bases=["A"],
)
# Returns "A" (high confidence)

# Or use Fitch parsimony directly
from ancify.ancestral import call_ancestral_base_parsimony
from ancify.parsimony import parse_newick

tree = parse_newick("(((bonobo,chimp),gorilla),macaque)")
base = call_ancestral_base_parsimony(tree, {
    "bonobo": "G", "chimp": "G", "gorilla": "A", "macaque": "A"
})
# Returns "A" (high confidence -- tree resolves the ambiguity)
```

---

## Evaluation (Optional)

Compare your calls against a reference ancestral sequence (e.g. Ensembl EPO) and/or VCF variant data:

```yaml
evaluation:
  reference_dir: ./ensembl_ancestor/
  reference_pattern: "homo_sapiens_ancestor_{chrom_id}.fa"
  vcf_dir: ./vcf/
  vcf_pattern: "ALL.chr{chrom_id}.vcf.gz"
```

Pattern placeholders: `{chrom}` = full name (e.g. `chr1`), `{chrom_id}` = without `chr` prefix (e.g. `1`).

### Human (hg38) Validation Results

The BCGM method (bonobo + chimp + gorilla + macaque) was validated against the Ensembl EPO 13-primate ancestral reference:

| Metric | chr1 | chr22 |
|--------|------|-------|
| Coverage (BCGM) | 79.1% | 74.6% |
| Coverage (Ensembl EPO) | 90.2% | 81.2% |
| Disagreement rate | 0.08% | 0.11% |
| Matches REF or ALT | 99.6% | 99.5% |

---

## Documentation

- **Online docs:** [ancify.readthedocs.io](https://ancify.readthedocs.io)

---

## Project Structure

```
ancify/
├── pyproject.toml              # package metadata
├── example_configs/
│   ├── hg38_bcgm.yaml         # human (worked example)
│   ├── mouse_example.yaml     # mouse (hypothetical)
│   ├── drosophila_example.yaml
│   └── brassica_rapa_example.yaml
├── docs/
├── tests/
└── ancify/                     # Python package
    ├── __init__.py
    ├── __main__.py
    ├── cli.py                  # command-line interface
    ├── config.py               # YAML config loading
    ├── utils.py                # FASTA I/O, majority vote
    ├── project.py              # Phase 1: coordinate projection
    ├── ancestral.py            # Phase 2: ancestral calling (dispatches by method)
    ├── parsimony.py            # Fitch algorithm & Newick parser
    ├── likelihood.py           # Felsenstein pruning, JC69/K80/HKY85/GTR
    ├── ml.py                   # LightGBM training & prediction
    ├── backend.py              # CPU/GPU vectorized voting
    └── evaluate.py             # Phase 3: evaluation
```

## License

MIT
