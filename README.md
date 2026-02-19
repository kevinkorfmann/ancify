# ancify

**Infer ancestral alleles for any species using outgroup alignments.**

ancify is a config-driven Python pipeline that determines the ancestral state at every position in a reference genome by comparing pairwise alignments from multiple outgroup species. It uses a two-tier inner/outer outgroup voting scheme with case-encoded confidence levels.

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
**Phase 2** infers the ancestral allele at every position using majority voting.
**Phase 3** (optional) evaluates calls against a reference and/or VCF variants.

### Confidence Encoding

| Character | Confidence | Meaning |
|-----------|-----------|---------|
| `ACGT` | High | Inner and outer outgroups agree |
| `acgt` | Low | Only one tier has data |
| `n` | Unresolved | Inner and outer disagree |
| `N` | Missing | No data from either tier |

---

## Installation

```bash
git clone https://github.com/kevinkorfmann/ancify.git
cd ancify
pip install .

# with evaluation dependencies (scikit-allel, matplotlib):
pip install '.[evaluate]'
```

**Requirements:** Python >= 3.8, PyYAML, NumPy.

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
```

Also works as: `python -m ancify run -c config.yaml`

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

See `example_configs/` for complete examples.

---

## How It Works

### The Algorithm

1. **Inner consensus**: majority vote among closely related outgroup species (e.g. bonobo, chimp, gorilla).
2. **Outer consensus**: majority vote among distantly related outgroup species (e.g. macaque).
3. **Compare**:
   - Agree &rarr; **high confidence** (uppercase)
   - One missing &rarr; **low confidence** (lowercase, use the available call)
   - Disagree &rarr; **unresolved** (`n`)
   - Both missing &rarr; **missing** (`N`)

This two-tier approach guards against incomplete lineage sorting and lineage-specific substitutions. The outer outgroup provides an independent evolutionary check on the inner consensus.

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

# Or call the core function directly
base = call_ancestral_base(
    inner_bases=["A", "A", "G"],
    outer_bases=["A"],
)
# Returns "A" (high confidence)
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

A comprehensive LaTeX manual with detailed algorithm descriptions, flowcharts, and worked examples is included in `docs/`. To compile:

```bash
cd docs && pdflatex manual.tex && pdflatex manual.tex
```

A pre-compiled PDF is included at `docs/manual.pdf`.

---

## Project Structure

```
ancify/
├── pyproject.toml              # package metadata
├── example_configs/
│   ├── hg38_bcgm.yaml         # human (worked example)
│   ├── mouse_example.yaml     # mouse (hypothetical)
│   └── drosophila_example.yaml # fruit fly (hypothetical)
├── docs/
│   ├── manual.tex              # comprehensive LaTeX manual
│   └── manual.pdf              # pre-compiled PDF
└── ancify/                     # Python package
    ├── __init__.py
    ├── __main__.py
    ├── cli.py                  # command-line interface
    ├── config.py               # YAML config loading
    ├── utils.py                # FASTA I/O, majority vote
    ├── project.py              # Phase 1: coordinate projection
    ├── ancestral.py            # Phase 2: ancestral calling
    └── evaluate.py             # Phase 3: evaluation
```

## License

MIT
