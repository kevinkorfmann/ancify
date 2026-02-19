# Evaluation

Phase 3 of the pipeline evaluates ancestral calls against external data. This
phase is **optional** and only runs if the `evaluation` block is present in the
configuration file.

## What gets evaluated

### Coverage statistics (always computed)

For each chromosome, ancify counts:

- **High-confidence positions** — uppercase `ACGT`.
- **Low-confidence positions** — lowercase `acgt`.
- **Missing positions** — `N` or `n`.
- **Proportions** — non-missing and high-confidence fractions.

### Reference comparison (optional)

If `reference_dir` and `reference_pattern` are configured, ancify compares its
calls against a reference ancestral sequence (e.g., Ensembl EPO). Metrics:

- **Agreement rate** — fraction of positions where both methods call and agree.
- **Disagreement rate** — fraction where both call but disagree.
- **Complementarity** — how often one method fills gaps left by the other.

### VCF comparison (optional)

If `vcf_dir` and `vcf_pattern` are configured, ancify checks how often its
ancestral call matches the REF or ALT allele at known variant sites. This
requires `scikit-allel`:

```bash
pip install 'ancify[evaluate]'
```

Metrics:

- **Proportion non-missing** — fraction of variant sites with an ancestral call.
- **Matches REF** — fraction matching the VCF reference allele.
- **Matches ALT** — fraction matching the VCF alternate allele.
- **Matches either** — should be >99% for a correct pipeline.

## Output format

Evaluation results are written as per-chromosome text files in
`<output_dir>/evaluation/<chrom>.evaluation.txt`:

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

## Validation results: Human hg38 (BCGM)

The BCGM method (bonobo + chimp + gorilla + macaque) was validated against the
Ensembl EPO 13-primate ancestral reference and 1000 Genomes Phase 3 VCFs:

| Metric | chr1 | chr22 | chrX |
|--------|------|-------|------|
| Coverage (BCGM) | 79.1% | 74.6% | 44.4% |
| Coverage (Ensembl EPO) | 90.2% | 81.2% | 44.1% |
| Disagreement rate | 0.08% | 0.11% | 3.51% |
| Matches REF/ALT | 99.6% | 99.5% | 99.3% |

Including low-confidence calls raises BCGM coverage from ~79% to ~90% on chr1,
with disagreement increasing from 0.08% to 0.57%.

## Configuring evaluation

```yaml
evaluation:
  reference_dir: ./ensembl_ancestor/
  reference_pattern: "homo_sapiens_ancestor_{chrom_id}.fa"
  vcf_dir: ./vcf/
  vcf_pattern: "ALL.chr{chrom_id}.vcf.gz"
```

Patterns support `{chrom}` (full name like `chr1`) and `{chrom_id}` (without
`chr` prefix, like `1`).

## Running evaluation standalone

```bash
ancify evaluate -c config.yaml
```

This only runs Phase 3; it assumes Phases 1 and 2 have already been completed
and the ancestral FASTA files exist in `output_dir`.
