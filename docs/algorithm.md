# Algorithm

This page describes the ancestral allele inference algorithm in detail.

## Overview

ancify infers the ancestral allele at every position in a focal species'
reference genome using a **two-tier outgroup voting scheme**:

1. **Inner outgroup consensus** — majority vote among closely related species.
2. **Outer outgroup consensus** — majority vote among distantly related species.
3. **Comparison** — the inner and outer consensuses are compared to produce a
   confidence-encoded ancestral call.

## Phase 1: Coordinate Projection

Before ancestral states can be called, each outgroup's pairwise alignment must
be converted into a sequence in the focal species' coordinate system.

### Input: Net AXT alignments

The pipeline reads **net AXT** pairwise alignment files from UCSC. These
represent best-in-genome one-to-one alignments between the focal and outgroup
species.

Each alignment block in an AXT file consists of four lines:

```
0 chr1 1 10 chrQ 500 509 + 1000
ACGTACGTNN          ← focal (target) sequence
ACGTACGTAA          ← outgroup (query) sequence
                    ← blank separator
```

### Projection logic

For each alignment block:

- Walk through the focal sequence character by character.
- If the focal character is a nucleotide (`A/T/C/G/N`), record the
  corresponding outgroup character at that focal position.
- If the focal character is a gap (`-`), skip — this represents an insertion
  in the outgroup that has no corresponding position in the focal genome.
- Positions not covered by any alignment block are filled with `N`.

The result is a FASTA file for each (outgroup, chromosome) pair, with the same
length as the focal chromosome.

## Phase 2: Ancestral State Inference

### Majority vote

For a set of bases from multiple species, the **majority vote** selects the
most frequent valid nucleotide (A, C, G, or T). Ties are broken alphabetically
(A > C > G > T). Bases that are `N`, `-`, or other non-nucleotide characters
are excluded from the vote.

If no base reaches the minimum frequency threshold (`min_inner_freq` or
`min_outer_freq`), the consensus is `N` (missing).

### Decision logic

At each genomic position, the algorithm:

1. Computes the **inner consensus** from all inner outgroup species.
2. Computes the **outer consensus** from all outer outgroup species.
3. Compares them:

| Inner | Outer | Result | Confidence |
|-------|-------|--------|------------|
| `A` | `A` | `A` | **High** (uppercase) |
| `A` | `N` | `a` | Low (lowercase) |
| `N` | `A` | `a` | Low (lowercase) |
| `A` | `T` | `n` | Unresolved (disagreement) |
| `N` | `N` | `N` | Missing (both tiers lack data) |

### Pseudocode

```python
def call_ancestral_base(inner_bases, outer_bases):
    inner = majority_vote(inner_bases)
    outer = majority_vote(outer_bases)

    if inner != "N" and inner == outer:
        return inner              # HIGH confidence (uppercase)
    if inner == "N" and outer != "N":
        return outer.lower()      # LOW confidence (lowercase)
    if inner != "N" and outer == "N":
        return inner.lower()      # LOW confidence (lowercase)
    if inner == "N" and outer == "N":
        return "N"                # Both missing
    return "n"                    # Disagreement
```

## Confidence encoding

The output FASTA encodes confidence via letter case:

| Character | Confidence | Condition |
|-----------|-----------|-----------|
| `ACGT` | High | Inner and outer outgroups agree |
| `acgt` | Low | Only one outgroup tier has data |
| `n` (lowercase) | Unresolved | Inner and outer disagree |
| `N` (uppercase) | Missing | Both tiers lack data |

This encoding is compatible with the Ensembl EPO ancestral sequence convention.

## The `min_inner_freq` parameter

This parameter controls how strict the inner majority vote is. With 3 inner
outgroup species:

| Value | Meaning |
|-------|---------|
| 1 | Any single species suffices (maximum coverage) |
| 2 | At least 2 of 3 must agree (balanced) |
| 3 | All 3 must agree (maximum stringency) |

The same logic applies to `min_outer_freq` for the outer outgroup tier.

## Biological rationale

### Why two tiers?

Using a single outgroup is vulnerable to:

- **Back mutations**: the outgroup lineage mutated away from the true ancestral state.
- **Incomplete lineage sorting (ILS)**: gene trees disagree with the species tree.
- **Alignment errors**: spurious alignments in repetitive regions.

The two-tier approach mitigates these risks:

- The **inner outgroup** (multiple closely related species) uses majority voting
  to overcome single-species errors.
- The **outer outgroup** (a more distant species) provides an independent
  evolutionary check. Agreement between the two tiers makes convergent
  misinference very unlikely.

### When the algorithm fails

The algorithm can still fail in rare cases:

- **Convergent substitution**: the same mutation occurred independently in both
  the focal and outer outgroup lineages.
- **ILS across the inner clade**: all inner species inherited the derived allele
  through ILS, and the outer outgroup independently carries the same derived allele.
- **Systematic alignment artifacts**: e.g., paralog confusion in segmental duplications.

These cases are flagged as either disagreement (`n`) or contribute to the
~0.1% error rate observed in practice.
