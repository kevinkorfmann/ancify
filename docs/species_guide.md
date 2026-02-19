# Adapting to Other Species

*ancify works with **any focal species** for which you have pairwise net AXT alignments. This guide teaches you how to think about outgroup selection and walks through examples from across the tree of life.*

---

## A framework for choosing outgroups

Before configuring ancify, you need to make one key decision: **which species are your inner outgroups and which are your outer outgroups?**

### Step 1: Draw the phylogeny

Sketch (or look up) the phylogenetic relationships around your focal species. You need at least three species total: your focal, one inner outgroup, and one outer outgroup.

```text
  Example for human:

                 6 Mya         9 Mya              25 Mya
                  ┌─── Bonobo
            ┌─────┤                          ┌ INNER tier
            │     └─── Chimp                 │ (closely related)
     ┌──────┤                                │
     │      └──────── Gorilla               ─┘
  ───┤
     │
     └──────────────────── Macaque          ── OUTER tier
                                               (distantly related)
```

### Step 2: Assign tiers

**Inner outgroups** should be:

- Closely related to the focal species (same genus or family)
- Diverged **more recently** than the outer outgroup
- Ideally 2 or more species (enables majority voting)

**Outer outgroups** should be:

- Clearly outside the inner clade
- Diverged at least 2-3x further than the inner outgroups
- Far enough that convergent mutations with the inner group are extremely rare

### Step 3: Check data availability

For each candidate outgroup, check if UCSC has a pairwise net AXT alignment to your focal assembly:

```text
https://hgdownload.soe.ucsc.edu/goldenPath/<focal_assembly>/
```

Look for directories named `vs<Outgroup>`. If the alignment exists, you are in business.

### Decision flowchart

```text
  Do you have ≥2 inner outgroups with AXT alignments?
    │
    ├─ YES → Great! Majority voting will be robust.
    │
    └─ NO → Do you have 1 inner + 1 outer?
              │
              ├─ YES → Still works. The outer outgroup provides
              │        the independent check. Consider adding more
              │        inner species if available.
              │
              └─ NO → You need at least 1 inner + 1 outer.
                      Check UCSC or generate your own alignments.
```

---

## Worked examples

### Human (hg38) — the gold standard

```yaml
focal_species: human
chromosome_lengths: chromoLens.txt

outgroups:
  inner:
    - name: bonobo      # ~6 Mya
      alignment: hg38.panPan3.net.axt.gz
    - name: chimp       # ~6 Mya
      alignment: hg38.panTro6.net.axt.gz
    - name: gorilla     # ~9 Mya
      alignment: hg38.gorGor6.net.axt.gz
  outer:
    - name: macaque     # ~25 Mya
      alignment: hg38.rheMac10.net.axt.gz

output_dir: ./human_ancestral
num_cpus: 24
```

**Why this works well:** Three inner outgroups provide redundancy. The inner-outer divergence ratio (~6-9 Mya vs. ~25 Mya) is large enough that convergent errors between tiers are negligible. Alignment coverage is excellent for all four species.

### Mouse (mm39) — minimal setup

```yaml
focal_species: mouse
chromosome_lengths: mm39.chromLens.txt

outgroups:
  inner:
    - name: rat            # ~12 Mya
      alignment: mm39.rn7.net.axt.gz
  outer:
    - name: rabbit         # ~90 Mya
      alignment: mm39.oryCun2.net.axt.gz

output_dir: ./mouse_ancestral
num_cpus: 8
```

**With only one inner species,** the majority vote is trivially that species' allele. The outer outgroup still provides the independent confirmation. To strengthen the inner tier, consider adding hamster or other rodents if alignments are available.

### *Drosophila melanogaster* (dm6) — non-chr naming

```yaml
focal_species: drosophila_melanogaster
chromosome_lengths: dm6.chromLens.txt
chromosomes: [2L, 2R, 3L, 3R, 4, X]

outgroups:
  inner:
    - name: simulans       # ~2.5 Mya
      alignment: dm6.droSim2.net.axt.gz
    - name: sechellia      # ~2.5 Mya
      alignment: dm6.droSec1.net.axt.gz
  outer:
    - name: yakuba         # ~6 Mya
      alignment: dm6.droYak3.net.axt.gz

output_dir: ./dmel_ancestral
num_cpus: 6
```

**Note:** The chromosome names (`2L`, `3R`, etc.) do not have a `chr` prefix — ancify handles any naming convention. The explicit `chromosomes` list excludes heterochromatic scaffolds.

### *Brassica rapa* (plant) — beyond animals

```yaml
focal_species: brassica_rapa
chromosome_lengths: braRap1.chromLens.txt

outgroups:
  inner:
    - name: brassica_oleracea   # close relative, same genus
      alignment: braRap1.braOleracea.net.axt.gz
  outer:
    - name: arabidopsis_thaliana  # ~20 Mya, same family Brassicaceae
      alignment: braRap1.araTha1.net.axt.gz

output_dir: ./brassica_rapa_ancestral
num_cpus: 4
```

**Plant genomes** often use chromosome naming like A01, A02. Omit the `chromosomes` key to process all entries in the lengths file. Plant genome alignments may be more fragmented due to whole-genome duplications — expect higher `N` rates than in mammals.

### Zebrafish (danRer11) — fish

```yaml
focal_species: zebrafish
chromosome_lengths: danRer11.chromLens.txt

outgroups:
  inner:
    - name: medaka
      alignment: danRer11.oryLat2.net.axt.gz
  outer:
    - name: fugu
      alignment: danRer11.fr3.net.axt.gz

output_dir: ./zebrafish_ancestral
num_cpus: 4
```

**Fish have high substitution rates** compared to mammals, so the divergence between inner and outer outgroups needs to be carefully considered. Medaka and fugu provide a reasonable tier separation for zebrafish.

---

## Getting the input data

### Net AXT alignments from UCSC

Download from:

```text
https://hgdownload.soe.ucsc.edu/goldenPath/<focal_assembly>/vs<Outgroup>/
```

Example for human vs. chimp:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro6/hg38.panTro6.net.axt.gz
```

### Chromosome lengths

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

---

## Tips for outgroup selection

### More inner species is better

Even 2 inner species is a major improvement over 1. With 1 inner species, a single alignment error or lineage-specific substitution produces a wrong call. With 2+, the majority vote provides robustness.

```text
  1 inner species:   accuracy ≈ alignment quality
  2 inner species:   accuracy ≈ max(alignment quality)
  3+ inner species:  accuracy ≈ consensus of multiple independent signals
```

### The outer outgroup must be clearly outside the inner clade

If the outer outgroup is too closely related to the inner species, ILS can affect both tiers *together*, producing false high-confidence calls:

```text
  BAD: outer is too close to inner

       ┌── Focal
    ┌──┤
    │  └── Inner 1             ILS can affect all three
  ──┤
    └──── Outer (barely        species → false agreement
          more distant)

  GOOD: outer is clearly distant

       ┌── Focal
    ┌──┤
    │  └── Inner 1
  ──┤
    │
    │
    └────────── Outer (deep   ILS between tiers is
                divergence)    negligible
```

### Stringency vs. coverage tradeoff

Increasing `min_inner_freq` requires more species to agree:

| Setting | Coverage | Accuracy |
|---------|----------|----------|
| `min_inner_freq=1` | Highest | Lower (but still >99%) |
| `min_inner_freq=2` (with 3 inner) | Moderate | High |
| `min_inner_freq=3` (with 3 inner) | Lowest | Highest |

For most demographic analyses, `min_inner_freq=1` is fine — the SFS shape is robust to rare misassignments.

### Sex chromosomes and other special cases

- **chrY** has very poor alignment coverage due to massive repetitive content. Expect >80% `N` (missing).
- **chrX** may have lower coverage than autosomes, especially in regions of reduced synteny.
- **Mitochondrial genomes** are not handled by standard UCSC pairwise alignments. For mtDNA polarization, consider using a multiple sequence alignment approach instead.
