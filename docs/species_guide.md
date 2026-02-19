# Adapting to Other Species

ancify works with **any focal species** for which you have pairwise net AXT
alignments to outgroup species.

## What you need

1. **A reference genome assembly** for your focal species.
2. **A chromosome-lengths file** — a tab-separated file with at least two
   columns (chromosome name, length in bp).
3. **Pairwise net AXT alignment files** for each outgroup species. Available
   from the [UCSC Genome Browser downloads](https://hgdownload.soe.ucsc.edu/downloads.html)
   for hundreds of species pairs.
4. **Phylogenetic knowledge** to decide which species are inner (closely related)
   and which are outer (distantly related).

## Choosing outgroups

### Inner outgroups (closely related)

- Should diverge from the focal species **more recently** than the outer outgroup.
- Using 2+ inner species enables majority voting, which is robust to single-species
  alignment gaps or lineage-specific substitutions.
- Ideal divergence: the same order or family.

### Outer outgroups (distantly related)

- Should be **far enough** that convergent substitutions with the inner outgroup
  are rare.
- A divergence of >2× the inner divergence is a good rule of thumb.
- Can also be multiple species with their own majority vote.

## Worked example: Human (hg38)

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

## Example: Mouse (mm39)

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

With only one inner species, the majority vote is trivially that species' allele.
The outer outgroup still provides the independent confirmation.

## Example: *Drosophila melanogaster* (dm6)

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

Note the non-`chr` chromosome naming (`2L`, `3R`, etc.) — ancify handles
any naming convention.

## Example: Zebrafish (danRer11)

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

## Getting the input data

### Net AXT alignments

Download from UCSC at:

```
https://hgdownload.soe.ucsc.edu/goldenPath/<focal_assembly>/vs<Outgroup>/
```

For example, human vs. chimp:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro6/hg38.panTro6.net.axt.gz
```

### Chromosome lengths

From a FASTA index:

```bash
samtools faidx reference.fa
cut -f1,2 reference.fa.fai > chromoLens.txt
```

From UCSC MySQL:

```bash
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A \
  -e "SELECT chrom, size FROM chromInfo" hg38 > chromoLens.txt
```

## Tips

- **More inner species is better** — it makes the majority vote more robust.
  Even 2 inner species is a significant improvement over 1.
- **The outer outgroup should be clearly outside the inner clade.** If the
  outer outgroup is too closely related to the inner species, ILS between the
  two tiers can cause false high-confidence calls.
- **Stringency vs. coverage tradeoff**: increase `min_inner_freq` for higher
  accuracy at the cost of more missing data.
- **chrY and other low-coverage chromosomes** will have high missingness rates
  regardless of method, due to repetitive content and poor alignment quality.
