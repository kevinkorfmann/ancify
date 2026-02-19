# Quickstart

*Go from zero to ancestral alleles in five minutes.*

**What you will learn:**

- How to install ancify
- How to configure it for your species
- How to run the pipeline
- How to read and use the output

**Prerequisites:** Python >= 3.8, a terminal, and pairwise net AXT alignment files for your outgroup species (see {doc}`species_guide` if you do not have these yet).

---

## Step 1: Install

```bash
git clone https://github.com/kevinkorfmann/ancify.git
cd ancify
pip install .
```

Verify it worked:

```bash
ancify --help
```

You should see a help message listing the available subcommands (`init`, `project`, `call`, `evaluate`, `run`).

## Step 2: Generate a config template

```bash
ancify init -o config.yaml
```

This creates a starter YAML file with all fields documented. Open it in your editor.

## Step 3: Edit the config

Fill in four things:

1. **`chromosome_lengths`** — a tab-separated file mapping chromosome names to their lengths.
2. **`outgroups.inner`** — one or more closely related species with their alignment files.
3. **`outgroups.outer`** — one or more distantly related species with their alignment files.
4. **`output_dir`** — where the ancestral FASTA files should go.

Here is a minimal human example:

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

:::{tip}
Not working with humans? See {doc}`species_guide` for ready-made configs for mouse, Drosophila, Brassica rapa, zebrafish, and guidance on setting up any other species.
:::

## Step 4: Run the pipeline

```bash
ancify run -c config.yaml
```

This runs all three phases in sequence:

| Phase | What it does | Time (human genome, 24 cores) |
|-------|-------------|-------------------------------|
| **1. Project** | Converts outgroup alignments into focal-species coordinates | 2-8 hours |
| **2. Call** | Infers the ancestral allele at every position | 5-15 minutes |
| **3. Evaluate** | Compares calls against a reference (optional) | 5-15 minutes |

:::{note}
Phase 1 is the slow step because it streams through large alignment files. Phases 2 and 3 are fast. If Phase 1 has already been run, you can re-run just Phase 2 with `ancify call -c config.yaml`.
:::

## Step 5: Understand the output

Your `output_dir` now contains one FASTA file per chromosome:

```text
ancestral_calls/
├── chr1.fa     ← 248 million characters, one per position
├── chr2.fa
├── ...
└── chrX.fa
```

Each position is a single character encoding both the ancestral call **and** its confidence:

| Character | Confidence | What it means |
|-----------|-----------|---------------|
| `A` `C` `G` `T` | **High** | Inner and outer outgroups agree — this is almost certainly the ancestral allele |
| `a` `c` `g` `t` | **Low** | Only one outgroup tier had data — probably correct but less certain |
| `n` | **Unresolved** | Inner and outer outgroups disagree — possible ILS or alignment artifact |
| `N` | **Missing** | No outgroup data at this position |

## Step 6: Use the output

### Look up a single position

```python
from ancify.utils import read_fasta

_, seq = read_fasta("ancestral_calls/chr1.fa")

position = 1_000_000           # 1-based genomic coordinate
allele = seq[position - 1]     # convert to 0-based index

print(f"Position {position}: ancestral = {allele}")
print(f"  High confidence: {allele.isupper() and allele in 'ACGT'}")
print(f"  Any call:        {allele.upper() in 'ACGT'}")
```

### Polarize variants from a VCF

```python
for variant in vcf:
    anc = seq[variant.POS - 1].upper()

    if anc not in "ACGT":
        continue  # no ancestral call at this position

    if anc == variant.REF:
        print(f"{variant.POS}: REF is ancestral, ALT is derived")
    elif anc == variant.ALT:
        print(f"{variant.POS}: ALT is ancestral, REF is derived")
    else:
        print(f"{variant.POS}: ancestral allele matches neither REF nor ALT")
```

### Filter by confidence

```python
if allele.isupper() and allele in "ACGT":
    # High confidence only — strictest, but most reliable
    use_for_analysis()
elif allele.upper() in "ACGT":
    # Include low confidence — more coverage, slightly less reliable
    use_for_analysis()
```

---

## What just happened?

Behind the scenes, ancify did three things:

```text
Phase 1: PROJECT
  For each outgroup species and each chromosome, it read the pairwise
  alignment file and projected the outgroup's bases onto the focal
  genome's coordinate system. This produced one FASTA file per
  (species, chromosome) pair, all with the same length as the focal
  chromosome.

Phase 2: CALL
  At every position, it collected the projected bases from all outgroup
  species, computed a majority vote among the inner outgroups and among
  the outer outgroups, compared the two, and wrote the ancestral call
  with confidence encoding.

Phase 3: EVALUATE (optional)
  If you provided a reference ancestral sequence and/or VCF files in
  the config, it compared the calls against these external data sources
  and wrote per-chromosome evaluation statistics.
```

For a deeper understanding of the algorithm, see {doc}`algorithm`. For the biology behind why this works, see {doc}`background`.

---

## Next steps

| I want to... | Go to... |
|--------------|----------|
| Understand the biology behind polarization | {doc}`background` |
| Walk through a complete worked example | {doc}`tutorials` |
| Configure ancify for a different species | {doc}`species_guide` |
| Understand every config field | {doc}`configuration` |
| Learn the algorithm in depth | {doc}`algorithm` |
| Interpret evaluation results | {doc}`evaluation` |
| Use ancify as a Python library | {doc}`api` |
| Look up a term | {doc}`glossary` |
