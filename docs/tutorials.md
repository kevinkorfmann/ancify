# Tutorials

*Hands-on, end-to-end walkthroughs that teach you both the tool and the science.*

---

## Tutorial 1: Polarizing the Human Genome

This is the original use case that motivated ancify. You will polarize all human chromosomes using four great-ape and Old World monkey outgroups.

### Background

Humans diverged from:
- **Bonobo & Chimpanzee** ~6 million years ago (Mya)
- **Gorilla** ~9 Mya
- **Macaque** ~25 Mya

These divergence times make bonobo, chimp, and gorilla excellent **inner outgroups** (close, high alignment coverage) and macaque a good **outer outgroup** (distant enough that shared derived alleles with the inner group are extremely rare).

```text
                 6 Mya        9 Mya             25 Mya
                  ┌─── Bonobo
            ┌─────┤                           INNER
            │     └─── Chimpanzee             outgroups
     ┌──────┤
     │      └──────── Gorilla
─────┤
     │
     └──────────────────── Macaque            OUTER outgroup
```

### Step 1: Download the data

You need four files from UCSC:

```bash
BASE=https://hgdownload.soe.ucsc.edu/goldenPath/hg38

wget $BASE/vsPanPan3/hg38.panPan3.net.axt.gz     # bonobo
wget $BASE/vsPanTro6/hg38.panTro6.net.axt.gz     # chimp
wget $BASE/vsGorGor6/hg38.gorGor6.net.axt.gz     # gorilla
wget $BASE/vsRheMac10/hg38.rheMac10.net.axt.gz   # macaque
```

And a chromosome lengths file:

```bash
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A \
  -e "SELECT chrom, size FROM chromInfo" hg38 > chromoLens.txt
```

### Step 2: Create the config

```yaml
focal_species: human
chromosome_lengths: chromoLens.txt

chromosomes:
  - chr1
  - chr2
  - chr3
  - chr4
  - chr5
  - chr6
  - chr7
  - chr8
  - chr9
  - chr10
  - chr11
  - chr12
  - chr13
  - chr14
  - chr15
  - chr16
  - chr17
  - chr18
  - chr19
  - chr20
  - chr21
  - chr22
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
output_dir: ./human_ancestor_bcgm
num_cpus: 24
```

:::{tip}
We deliberately omit chrY here. The Y chromosome has very poor alignment coverage due to its repetitive content, so the ancestral calls would be mostly `N` (missing). Including it is harmless but not very informative.
:::

### Step 3: Run the pipeline

```bash
ancify run -c hg38_bcgm.yaml
```

**Expected runtime** on a machine with 24 cores and sufficient RAM (~32 GB):

| Phase | Time |
|-------|------|
| Phase 1 (projection) | 2-8 hours |
| Phase 2 (calling) | 5-15 minutes |

:::{note}
Phase 1 is I/O-bound — it streams through ~6 GB of compressed alignment data. If you have an SSD, it will be faster. If you want to monitor progress, use verbose mode: `ancify -v run -c hg38_bcgm.yaml`.
:::

### Step 4: Inspect the results

```python
from ancify.utils import read_fasta

_, seq = read_fasta("human_ancestor_bcgm/chr22.fa")

total = len(seq)
high = sum(1 for c in seq if c in "ACGT")
low = sum(1 for c in seq if c in "acgt")
unresolved = sum(1 for c in seq if c == "n")
missing = sum(1 for c in seq if c == "N")

print(f"Chr22 length:     {total:>12,}")
print(f"High confidence:  {high:>12,}  ({high/total:.1%})")
print(f"Low confidence:   {low:>12,}  ({low/total:.1%})")
print(f"Unresolved:       {unresolved:>12,}  ({unresolved/total:.1%})")
print(f"Missing:          {missing:>12,}  ({missing/total:.1%})")
```

Typical output for chr22:

```text
Chr22 length:        50,818,468
High confidence:     37,910,573  (74.6%)
Low confidence:       8,021,119  (15.8%)
Unresolved:              34,562  (0.1%)
Missing:              4,852,214  (9.5%)
```

### What do these numbers mean?

- **74.6% high confidence** means that at nearly three-quarters of all positions, all four outgroups provided data and the inner and outer tiers agreed. These are your most reliable calls.
- **15.8% low confidence** means only one tier had data (typically the inner outgroups had coverage but macaque did not, or vice versa). These are still useful but less certain.
- **0.1% unresolved** means the inner and outer tiers disagreed. This could indicate incomplete lineage sorting, a back-mutation in one lineage, or an alignment artifact. It is a tiny fraction, which is reassuring.
- **9.5% missing** means no outgroup had alignable sequence at these positions — typically repetitive regions, centromeres, or assembly gaps.

### Step 5: Add evaluation (optional)

To validate against the Ensembl EPO ancestral reference:

```yaml
evaluation:
  reference_dir: ./ensembl_ancestor/
  reference_pattern: "homo_sapiens_ancestor_{chrom_id}.fa"
  vcf_dir: ./1kg_vcf/
  vcf_pattern: "ALL.chr{chrom_id}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
```

Run evaluation alone (if Phases 1-2 already completed):

```bash
ancify evaluate -c hg38_bcgm.yaml
```

See {doc}`evaluation` for how to interpret the results.

---

## Tutorial 2: Your First Non-Human Species

This tutorial walks through the thought process of setting up ancify for a new species, using mouse as an example.

### Think phylogenetically

Before touching the keyboard, draw (or look up) the phylogeny around your focal species:

```text
                12 Mya                          90 Mya
                  ┌──── Mouse (focal)
     ┌────────────┤
     │            └──── Rat (inner outgroup)
─────┤
     │
     └──────────────────────── Rabbit (outer outgroup)
```

**Ask yourself:**

1. **Which species are closely related?** These become your inner outgroups. Rat diverged from mouse ~12 Mya — close enough for good alignment coverage.
2. **Which species are clearly outside that clade?** These become your outer outgroup. Rabbit diverged ~90 Mya — far enough that convergent errors between the tiers are negligible.
3. **Are pairwise alignments available?** Check UCSC at `https://hgdownload.soe.ucsc.edu/goldenPath/mm39/`.

### Get the data

```bash
# Alignment files
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsRn7/mm39.rn7.net.axt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsOryCun2/mm39.oryCun2.net.axt.gz

# Chromosome lengths
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A \
  -e "SELECT chrom, size FROM chromInfo" mm39 > mm39.chromLens.txt
```

### Create the config

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
num_cpus: 8
```

:::{important}
With only one inner species, the "majority vote" is trivially that species' allele. The outer outgroup still provides the independent check. For more robust calls, consider adding more inner species (e.g. hamster if a UCSC alignment is available).
:::

### Run and inspect

```bash
ancify run -c mouse_config.yaml
```

The output follows the same confidence encoding as any other species. See {doc}`species_guide` for more examples.

---

## Tutorial 3: Interpreting Disagreements

Not all positions yield a clean ancestral call. This tutorial teaches you what the edge cases mean and how to handle them.

### Case 1: High confidence (`A`)

```text
  Bonobo:   A
  Chimp:    A     → Inner majority: A
  Gorilla:  A
  Macaque:  A     → Outer: A

  Inner == Outer → "A" (high confidence, uppercase)
```

**Interpretation:** Four independent species all carry A. Under parsimony, A is almost certainly the ancestral allele. The probability of all four independently mutating to A from something else is vanishingly small.

### Case 2: Low confidence, inner only (`a`)

```text
  Bonobo:   A
  Chimp:    A     → Inner majority: A
  Gorilla:  A
  Macaque:  N     → Outer: N (no alignment at this position)

  Inner has data, Outer missing → "a" (low confidence, lowercase)
```

**Interpretation:** The inner outgroups agree on A, but we have no independent confirmation from the outer outgroup. This is probably correct — three species agreeing is strong evidence — but we cannot rule out ILS affecting the entire inner clade.

### Case 3: Low confidence, outer only (`t`)

```text
  Bonobo:   N
  Chimp:    N     → Inner majority: N (all missing)
  Gorilla:  N
  Macaque:  T     → Outer: T

  Inner missing, Outer has data → "t" (low confidence, lowercase)
```

**Interpretation:** Only the distant outgroup has data. This typically happens in regions where the focal species has diverged enough from the inner outgroups to lose alignment, but the more distant outgroup still aligns (perhaps due to a conserved element). Use with caution.

### Case 4: Disagreement (`n`)

```text
  Bonobo:   A
  Chimp:    A     → Inner majority: A
  Gorilla:  G
  Macaque:  G     → Outer: G

  Inner says A, Outer says G → "n" (unresolved)
```

**Interpretation:** This is the interesting case. Possible explanations:

1. **Incomplete lineage sorting (ILS):** The gene tree at this position disagrees with the species tree. The inner species inherited one allele from the ancestral population, macaque inherited the other.
2. **Convergent mutation:** A mutation happened independently in the inner lineage and the outer lineage arrived at different alleles.
3. **Alignment artifact:** One of the alignments is incorrect at this position.

ancify flags these as `n` rather than guessing. For most analyses, excluding these positions (typically <0.2% of the genome) is the safest approach.

### Case 5: Missing (`N`)

```text
  Bonobo:   N
  Chimp:    N     → Inner majority: N
  Gorilla:  N
  Macaque:  N     → Outer: N

  Both tiers missing → "N"
```

**Interpretation:** No outgroup has alignable sequence here. Common in centromeres, telomeres, segmental duplications, and transposable element insertions. Nothing can be inferred at these positions.

---

## Tutorial 4: Getting Ancestral FASTA Files

After running the pipeline you have one FASTA file per chromosome in your
`output_dir`. This tutorial shows how to verify, inspect, and use them.

### Locate the output

The output directory mirrors the `output_dir` field in your config:

```bash
ls -lh ancestral_calls/
```

```text
ancestral_calls/
├── chr1.fa      248 MB
├── chr2.fa      242 MB
├── ...
├── chr22.fa      51 MB
└── chrX.fa      156 MB
```

Each file is a plain-text FASTA with a single record whose sequence is
exactly as long as the corresponding chromosome.

### Quick sanity check

Verify that the file length matches the chromosome length:

```bash
# Print the sequence length (excluding the header line)
awk 'NR>1{len+=length($0)} END{print len}' ancestral_calls/chr22.fa
# Expected: 50818468 (hg38 chr22 length)
```

Or from Python:

```python
from ancify.utils import read_fasta

_, seq = read_fasta("ancestral_calls/chr22.fa")
print(f"chr22 length: {len(seq):,}")
```

### Read a single position

Positions are 0-indexed in the sequence string (matching Python conventions)
but 1-indexed in genomic coordinates. To look up a 1-based genomic position:

```python
from ancify.utils import read_fasta

_, seq = read_fasta("ancestral_calls/chr22.fa")

pos = 16_050_075                # 1-based genomic coordinate
allele = seq[pos - 1]           # convert to 0-based

if allele.upper() in "ACGT":
    confidence = "high" if allele.isupper() else "low"
    print(f"chr22:{pos}  ancestral = {allele.upper()}  ({confidence} confidence)")
else:
    print(f"chr22:{pos}  no ancestral call ({allele})")
```

### Load all chromosomes into a dictionary

For genome-wide analyses it is convenient to load everything once:

```python
from pathlib import Path
from ancify.utils import read_fasta

anc_dir = Path("ancestral_calls")
ancestral = {}

for fa in sorted(anc_dir.glob("chr*.fa")):
    chrom = fa.stem                       # "chr1", "chr22", ...
    _, ancestral[chrom] = read_fasta(fa)
    print(f"  loaded {chrom} ({len(ancestral[chrom]):,} bp)")

print(f"\nLoaded {len(ancestral)} chromosomes")
```

### Confidence encoding recap

| Character | Meaning |
|-----------|---------|
| `A` `C` `G` `T` | High confidence — inner and outer outgroups agree |
| `a` `c` `g` `t` | Low confidence — only one tier had data |
| `n` | Unresolved — the two tiers disagree |
| `N` | Missing — no outgroup data |

Filter by confidence depending on your tolerance for noise:

```python
# Strict: high-confidence only
if allele in "ACGT":
    use(allele)

# Lenient: any call (high or low)
if allele.upper() in "ACGT":
    use(allele.upper())
```

---

## Tutorial 5: Polarizing VCF Variants

*Annotate every variant in a VCF as ancestral or derived.*

Polarization means labelling which allele (REF or ALT) was present in the
common ancestor. This is essential for computing the **unfolded site
frequency spectrum**, McDonald–Kreitman tests, and any analysis that
distinguishes ancestral from derived mutations.

### Prerequisites

- Ancestral FASTA files from Tutorial 4 (or the quickstart)
- A VCF file for the same reference genome (e.g. 1000 Genomes on GRCh38)
- `cyvcf2` (recommended) or `pysam` for reading VCFs — install with
  `pip install cyvcf2`

### Concept

For each biallelic SNP in the VCF:

1. Look up the ancestral allele at that position from the FASTA.
2. Compare it to REF and ALT.
3. If the ancestral allele matches REF → ALT is derived.
4. If the ancestral allele matches ALT → REF is derived (the reference
   genome carries the derived allele at this site).
5. If it matches neither → the variant cannot be polarized.

### Minimal example

```python
from cyvcf2 import VCF
from ancify.utils import read_fasta

_, seq = read_fasta("ancestral_calls/chr22.fa")

for v in VCF("chr22.vcf.gz"):
    if len(v.ALT) != 1 or len(v.REF) != 1 or len(v.ALT[0]) != 1:
        continue  # skip multiallelic and indels

    anc = seq[v.POS - 1].upper()
    if anc not in "ACGT":
        continue  # no ancestral call

    if anc == v.REF:
        print(f"{v.CHROM}:{v.POS}  REF={v.REF} is ancestral, ALT={v.ALT[0]} is derived")
    elif anc == v.ALT[0]:
        print(f"{v.CHROM}:{v.POS}  ALT={v.ALT[0]} is ancestral, REF={v.REF} is derived")
    else:
        print(f"{v.CHROM}:{v.POS}  ancestral {anc} matches neither REF nor ALT")
```

### Complete script: write an annotated VCF

The script below reads a VCF, adds an `AA` INFO field (ancestral allele),
and writes a new VCF with polarization annotations:

```python
#!/usr/bin/env python3
"""Polarize a VCF using ancify ancestral FASTA files.

Usage:
    python polarize_vcf.py \
        --anc-dir ancestral_calls/ \
        --vcf input.vcf.gz \
        --out polarized.vcf.gz
"""

import argparse
from pathlib import Path

from cyvcf2 import VCF, Writer
from ancify.utils import read_fasta


def load_ancestral(anc_dir):
    """Load all ancestral FASTA files into a dict keyed by chromosome name."""
    seqs = {}
    for fa in Path(anc_dir).glob("*.fa"):
        chrom = fa.stem
        _, seqs[chrom] = read_fasta(fa)
    return seqs


def polarize(vcf_path, anc_dir, out_path, high_confidence_only=False):
    vcf = VCF(str(vcf_path))
    vcf.add_info_to_header({
        "ID": "AA",
        "Number": "1",
        "Type": "String",
        "Description": "Ancestral allele inferred by ancify",
    })
    vcf.add_info_to_header({
        "ID": "DAF",
        "Number": "1",
        "Type": "String",
        "Description": "Derived allele: REF or ALT (or . if ambiguous)",
    })

    w = Writer(str(out_path), vcf)
    ancestral = load_ancestral(anc_dir)

    counts = {"ref_ancestral": 0, "alt_ancestral": 0, "nomatch": 0, "skipped": 0}

    for v in vcf:
        seq = ancestral.get(v.CHROM)
        if seq is None or v.POS - 1 >= len(seq):
            counts["skipped"] += 1
            w.write_record(v)
            continue

        allele = seq[v.POS - 1]
        if high_confidence_only and not allele.isupper():
            counts["skipped"] += 1
            w.write_record(v)
            continue

        anc = allele.upper()
        if anc not in "ACGT":
            counts["skipped"] += 1
            w.write_record(v)
            continue

        v.INFO["AA"] = anc

        if len(v.ALT) == 1 and len(v.REF) == 1 and len(v.ALT[0]) == 1:
            if anc == v.REF:
                v.INFO["DAF"] = "ALT"
                counts["ref_ancestral"] += 1
            elif anc == v.ALT[0]:
                v.INFO["DAF"] = "REF"
                counts["alt_ancestral"] += 1
            else:
                v.INFO["DAF"] = "."
                counts["nomatch"] += 1

        w.write_record(v)

    w.close()
    vcf.close()

    total = sum(counts.values())
    print(f"Polarized {total:,} variants:")
    print(f"  REF is ancestral:  {counts['ref_ancestral']:>10,}")
    print(f"  ALT is ancestral:  {counts['alt_ancestral']:>10,}")
    print(f"  No match:          {counts['nomatch']:>10,}")
    print(f"  Skipped:           {counts['skipped']:>10,}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Polarize VCF with ancestral alleles")
    parser.add_argument("--anc-dir", required=True, help="Directory with ancestral chr*.fa files")
    parser.add_argument("--vcf", required=True, help="Input VCF (optionally gzipped)")
    parser.add_argument("--out", required=True, help="Output VCF path")
    parser.add_argument("--high-only", action="store_true",
                        help="Only use high-confidence calls (uppercase)")
    args = parser.parse_args()
    polarize(args.vcf, args.anc_dir, args.out, high_confidence_only=args.high_only)
```

Save this as `polarize_vcf.py` and run:

```bash
python polarize_vcf.py \
    --anc-dir ancestral_calls/ \
    --vcf ALL.chr22.vcf.gz \
    --out chr22_polarized.vcf.gz
```

### Using scikit-allel instead of cyvcf2

If you already have `scikit-allel` installed (it ships with
`pip install 'ancify[evaluate]'`), you can polarize without cyvcf2:

```python
import allel
import numpy as np
from ancify.utils import read_fasta

_, seq = read_fasta("ancestral_calls/chr22.fa")

vcf = allel.read_vcf(
    "chr22.vcf.gz",
    fields=["variants/POS", "variants/REF", "variants/ALT"],
)

pos = vcf["variants/POS"]
ref = vcf["variants/REF"]
alt = vcf["variants/ALT"][:, 0]

anc = np.array(list(seq), dtype="U1")[pos - 1]
anc_upper = np.char.upper(anc)

valid = np.isin(anc_upper, list("ACGT"))
ref_is_anc = (anc_upper == ref) & valid
alt_is_anc = (anc_upper == alt) & valid

print(f"Total SNPs:          {len(pos):>10,}")
print(f"With ancestral call: {valid.sum():>10,}")
print(f"REF is ancestral:    {ref_is_anc.sum():>10,}")
print(f"ALT is ancestral:    {alt_is_anc.sum():>10,}")
```

### Computing the unfolded site frequency spectrum

Once variants are polarized, computing the unfolded SFS is straightforward.
The derived allele frequency at each site is the frequency of whichever
allele is *not* the ancestral one:

```python
import allel
import numpy as np
from ancify.utils import read_fasta

_, seq = read_fasta("ancestral_calls/chr22.fa")

vcf = allel.read_vcf(
    "chr22.vcf.gz",
    fields=["variants/POS", "variants/REF", "variants/ALT", "calldata/GT"],
)

pos = vcf["variants/POS"]
ref = vcf["variants/REF"]
alt = vcf["variants/ALT"][:, 0]
gt = allel.GenotypeArray(vcf["calldata/GT"])

anc = np.array(list(seq), dtype="U1")[pos - 1]
anc_upper = np.char.upper(anc)

ref_is_anc = anc_upper == ref
alt_is_anc = anc_upper == alt
valid = ref_is_anc | alt_is_anc

ac = gt[valid].count_alleles()
n_samples = gt.shape[1]
n_chroms = 2 * n_samples

dac = np.where(
    ref_is_anc[valid],
    ac[:, 1],          # ALT count = derived count
    ac[:, 0],          # REF count = derived count (ALT is ancestral)
)

sfs = np.bincount(dac, minlength=n_chroms + 1)

print("Unfolded SFS (first 10 bins):")
for i in range(min(10, len(sfs))):
    print(f"  DAC={i}: {sfs[i]:,} sites")
```

### Tips

- **Use high-confidence calls for SFS analyses.** Low-confidence calls can
  introduce polarization errors that inflate the high-frequency derived
  class. Pass `--high-only` to the script above, or filter in code with
  `allele.isupper()`.
- **Exclude `n` (unresolved) sites.** These are positions where inner and
  outer outgroups disagree — possibly due to incomplete lineage sorting.
  They are already excluded by the `anc not in "ACGT"` check above.
- **Watch for REF bias.** In most VCFs, REF is the reference genome allele,
  which is often the major allele. Expect `REF is ancestral` to be the most
  common outcome (~90% for human common variants).

---

## Next steps

| I want to... | Go to... |
|--------------|----------|
| Set up a new species from scratch | {doc}`species_guide` |
| Understand every config option | {doc}`configuration` |
| Dive deep into the algorithm | {doc}`algorithm` |
| Learn what the evaluation numbers mean | {doc}`evaluation` |
| Fix something that is not working | {doc}`faq` |
