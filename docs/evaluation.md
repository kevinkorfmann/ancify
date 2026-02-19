# Evaluation

*Phase 3 of the pipeline evaluates ancestral calls against external data. This page explains what gets measured, how to interpret the results, and what the numbers mean for your analysis.*

---

## Overview

Evaluation is **optional** and only runs if the `evaluation` block is present in your configuration file. It answers two questions:

1. **How do my calls compare to a known reference?** (e.g., Ensembl EPO ancestral sequence)
2. **Do my calls make sense at known variant sites?** (e.g., 1000 Genomes VCF)

---

## What gets measured

### Coverage statistics (always computed)

For each chromosome, ancify counts:

| Metric | What it measures |
|--------|-----------------|
| High-confidence positions | Uppercase `ACGT` — both tiers agree |
| Low-confidence positions | Lowercase `acgt` — one tier has data |
| Unresolved positions | `n` — tiers disagree |
| Missing positions | `N` — no data |

### Reference comparison (optional)

If `reference_dir` and `reference_pattern` are configured, ancify compares its calls position-by-position against a reference ancestral sequence.

| Metric | What it measures |
|--------|-----------------|
| Agreement rate | Fraction of positions where both methods call and agree |
| Disagreement rate | Fraction where both call but disagree |
| Complementarity | How often one method fills gaps left by the other |

### VCF comparison (optional)

If `vcf_dir` and `vcf_pattern` are configured, ancify checks how often its ancestral call matches the REF or ALT allele at known variant sites.

| Metric | What it measures |
|--------|-----------------|
| Proportion non-missing | Fraction of variant sites with an ancestral call |
| Matches REF | Fraction matching the VCF reference allele |
| Matches ALT | Fraction matching the VCF alternate allele |
| Matches either | Should be >99% for a correct pipeline |

---

## Output format

Results are written as per-chromosome text files in `<output_dir>/evaluation/`:

```ini
[coverage]
  total_positions: 248956422
  high_confidence: 197145832
  low_confidence: 25203651
  missing: 26606939
  prop_nonmissing: 0.8931
  prop_high_confidence: 0.7918

[reference_comparison]
  test_nonmissing: 0.7912
  ref_nonmissing: 0.9016
  both_nonmissing: 5234112
  agreement_rate: 0.9992
  disagreement_rate: 0.0008

[vcf_comparison]
  num_variants: 6102435
  prop_nonmissing: 0.8965
  matches_ref: 0.9453
  matches_alt: 0.0494
  matches_either: 0.9948
```

---

## How to interpret the results

### Coverage: "How much of my genome has an ancestral call?"

```text
  ┌─────────────────────────────────────────────────────────┐
  │                    Total genome                          │
  │                                                         │
  │  ████████████████████████████  ░░░░░░░░  ▓▓▓▓  ████    │
  │  HIGH confidence (75%)         LOW (16%)  n(0.1) N(9%)  │
  │                                                         │
  └─────────────────────────────────────────────────────────┘
```

**Typical values for human (autosomes):**

- `prop_nonmissing` ≈ 0.89-0.91 (89-91% of positions have some call)
- `prop_high_confidence` ≈ 0.75-0.80 (75-80% are high confidence)

**If your values are much lower:**

- Check that all alignment files were downloaded completely (`gzip -t file.axt.gz`)
- Verify chromosome names match between your lengths file and alignment files
- More distant outgroups will naturally have lower coverage

### Agreement rate: "Do I agree with the gold standard?"

```text
  agreement_rate:    0.9992    ← 99.92% of compared positions match
  disagreement_rate: 0.0008   ← 0.08% differ
```

**What does 0.08% disagreement mean?**

On a chromosome with 200 million callable positions, 0.08% is ~160,000 positions. This sounds like a lot, but in the context of the SFS:

- Most disagreements are at positions with rare alleles or in regions of poor alignment quality.
- The SFS is a *summary statistic* — individual site errors are diluted across frequency bins.
- Demographic inference tools are robust to error rates well below 1%.

**If your disagreement rate is much higher:**

- Check alignment quality for the most divergent outgroup
- Consider increasing `min_inner_freq` for stricter consensus
- Examine whether the disagreements are clustered in specific genomic regions (e.g., segmental duplications)

### VCF comparison: "Does my ancestral allele match known variants?"

```text
  matches_ref:    0.9453    ← 94.5% of variant sites: ancestral = REF
  matches_alt:    0.0494    ←  4.9% of variant sites: ancestral = ALT
  matches_either: 0.9948    ← 99.5% match one or the other
```

**Interpreting these numbers:**

- **`matches_ref` ≈ 95%** is expected. For most variants, the reference genome carries the ancestral allele (by construction — reference assemblies are built from common alleles).
- **`matches_alt` ≈ 5%** means 5% of variant sites have a derived reference allele. These are the sites where polarization matters most — it tells you to flip the frequency.
- **`matches_either` > 99%** is the key quality metric. This should always be very high. Values below 99% suggest a problem with the ancestral calls or the VCF data.

**What about the ~0.5% that match neither?**

These are positions where the ancestral call is a base that is neither the VCF REF nor ALT. Causes:

- Triallelic sites (the VCF REF and ALT are both derived, the ancestral allele is a third base)
- VCF errors or multi-allelic sites recorded as biallelic
- Ancestral calling errors (rare)

---

## Validation results: Human hg38 (BCGM)

The BCGM method (bonobo + chimp + gorilla + macaque) was validated against the Ensembl EPO 13-primate ancestral reference and 1000 Genomes Phase 3 VCFs:

| Metric | chr1 | chr22 | chrX |
|--------|------|-------|------|
| Coverage (BCGM) | 79.1% | 74.6% | 44.4% |
| Coverage (Ensembl EPO) | 90.2% | 81.2% | 44.1% |
| Disagreement rate | 0.08% | 0.11% | 3.51% |
| Matches REF/ALT | 99.6% | 99.5% | 99.3% |

**Key observations:**

- **Autosomes** show excellent agreement (<0.12% disagreement).
- **chrX** has higher disagreement (3.51%), reflecting male hemizygosity and reduced effective population size, which increases ILS rates.
- **Coverage** is lower than Ensembl EPO because ancify uses 4 species vs. EPO's 13. Including low-confidence calls raises BCGM coverage from ~79% to ~90% on chr1.

---

## Configuring evaluation

```yaml
evaluation:
  reference_dir: ./ensembl_ancestor/
  reference_pattern: "homo_sapiens_ancestor_{chrom_id}.fa"
  vcf_dir: ./vcf/
  vcf_pattern: "ALL.chr{chrom_id}.vcf.gz"
```

You can include either, both, or neither of `reference_dir` and `vcf_dir`:

- **Reference only:** Compare against a known ancestral sequence.
- **VCF only:** Check that ancestral calls match known variants.
- **Both:** Full validation suite.
- **Neither:** Evaluation block present but empty — only coverage statistics are computed.

### Pattern placeholders

| Placeholder | Example (`chr1`) | Example (`2L`) |
|-------------|-----------------|----------------|
| `{chrom}` | `chr1` | `2L` |
| `{chrom_id}` | `1` | `2L` |

---

## Running evaluation standalone

If Phases 1 and 2 have already been completed:

```bash
ancify evaluate -c config.yaml
```

This only runs Phase 3, reading the existing ancestral FASTA files from `output_dir`.

---

## When you do not need evaluation

Evaluation requires external reference data, which may not exist for non-model organisms. In that case:

- Omit the `evaluation` block from your config entirely.
- Phase 3 will be skipped.
- You can still assess quality by inspecting coverage statistics manually (count uppercase, lowercase, `n`, and `N` in the output FASTAs).

The {doc}`tutorials` page includes Python code for computing these statistics from the output files.
