# Algorithm

*A deep dive into how ancify infers ancestral alleles. This page teaches you the algorithm from first principles, with worked examples at every step.*

---

## Overview

ancify infers the ancestral allele at every position in a focal species' reference genome using a **two-tier outgroup voting scheme**:

```text
  INPUT                           PROCESSING                        OUTPUT
  ─────                           ──────────                        ──────

  Net AXT files        ┌──────────────────────┐
  (one per outgroup    │  Phase 1: PROJECT    │    Projected FASTAs
   species)       ────▶│  Map outgroup bases  │───▶ (one per species
                       │  to focal coords     │     per chromosome)
                       └──────────────────────┘
                                                         │
                                                         ▼
                       ┌──────────────────────┐
                       │  Phase 2: CALL       │    Ancestral FASTAs
                  ────▶│  Two-tier voting     │───▶ (one per chromosome,
                       │  + confidence encode │     confidence-encoded)
                       └──────────────────────┘
```

---

## Phase 1: Coordinate Projection

Before ancestral states can be called, each outgroup's pairwise alignment must be converted into a sequence in the focal species' coordinate system. This is the most computationally expensive step.

### Input: Net AXT alignments

The pipeline reads **net AXT** pairwise alignment files from UCSC. These represent best-in-genome one-to-one alignments between the focal and outgroup species.

Each alignment block in an AXT file consists of four lines:

```text
  Line 1: header    0 chr1 100 110 chrQ 500 510 + 1000
  Line 2: target    ACGTNNACGT       ← focal species' sequence
  Line 3: query     ACGTAAACGT       ← outgroup species' sequence
  Line 4:                            ← blank separator
```

**Why net alignments?** Raw pairwise alignments can have overlapping blocks (e.g. from segmental duplications). The "netting" algorithm resolves these overlaps by keeping only the highest-scoring chain at each position, producing clean one-to-one mappings.

### How projection works

For each alignment block, the algorithm walks through the focal sequence character by character:

```text
  Focal (target):   A  C  G  -  T  N  N  A  C  G  T
  Outgroup (query): A  C  G  T  T  A  A  A  C  G  T
                    ↓  ↓  ↓     ↓  ↓  ↓  ↓  ↓  ↓  ↓
  Projected:        A  C  G     T  A  A  A  C  G  T

  Rules:
    - Focal is a base → record the outgroup's base at that focal position
    - Focal is a gap (-) → skip (insertion in outgroup, no focal coordinate)
    - Position not covered by any block → fill with N
```

The result is a FASTA file for each (outgroup, chromosome) pair, with the **same length** as the focal chromosome. Every position either contains the outgroup's aligned base or `N` (no alignment).

### Worked example

Suppose the focal chromosome is 20 bases long and we have two alignment blocks:

```text
  Focal chromosome:    positions 1-20 (20 bases)

  Block A (positions 3-8):
    Focal:    A C G T A C
    Outgroup: A C A T A C

  Block B (positions 15-18):
    Focal:    T G A C
    Outgroup: T G A T

  Projected sequence:
    Position:  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20
    Result:    N  N  A  C  A  T  A  C  N  N  N  N  N  N  T  G  A  T  N  N
                     └── Block A ──┘                     └─ Block B ┘
```

Positions 1-2, 9-14, and 19-20 have no alignment coverage and are filled with `N`.

---

## Phase 2: Ancestral State Inference

This is the core algorithm. At every genomic position, it collects the projected bases from all outgroups, computes votes, and produces a confidence-encoded ancestral call.

### Step 1: Majority vote

For a set of bases from multiple species, the **majority vote** selects the most frequent valid nucleotide (A, C, G, or T):

```text
  Example: inner outgroups at position chr1:1000

    Bonobo:   A
    Chimp:    A
    Gorilla:  G

    Tally:  A=2, G=1
    Majority vote → A (frequency 2)
```

**Tie-breaking:** Ties are broken alphabetically (A > C > G > T). This is arbitrary but deterministic.

**Filtering:** Bases that are `N`, `-`, or other non-nucleotide characters are excluded from the vote. If no base reaches the minimum frequency threshold (`min_inner_freq` or `min_outer_freq`), the consensus is `N` (missing).

```text
  Example with min_inner_freq=2:

    Bonobo:   A
    Chimp:    N       ← excluded (missing)
    Gorilla:  G

    Tally:  A=1, G=1
    Neither reaches threshold of 2 → consensus is N
```

### Step 2: Compare inner and outer

Once the inner consensus and outer consensus are computed, the algorithm compares them:

```text
  ┌─────────┬─────────┬─────────────┬─────────────────────────────────┐
  │  Inner  │  Outer  │   Output    │  Interpretation                 │
  ├─────────┼─────────┼─────────────┼─────────────────────────────────┤
  │    A    │    A    │   A (HIGH)  │  Both tiers agree → confident   │
  │    A    │    N    │   a (low)   │  Only inner has data            │
  │    N    │    A    │   a (low)   │  Only outer has data            │
  │    A    │    T    │   n (none)  │  Tiers disagree → unresolved    │
  │    N    │    N    │   N (none)  │  No data at all                 │
  └─────────┴─────────┴─────────────┴─────────────────────────────────┘
```

### The complete algorithm as pseudocode

```python
def call_ancestral_base(inner_bases, outer_bases,
                        min_inner_freq=1, min_outer_freq=1):
    inner = majority_vote(inner_bases, min_freq=min_inner_freq)
    outer = majority_vote(outer_bases, min_freq=min_outer_freq)

    if inner != "N" and inner == outer:
        return inner              # HIGH confidence (uppercase)
    if inner == "N" and outer != "N":
        return outer.lower()      # LOW confidence (lowercase)
    if inner != "N" and outer == "N":
        return inner.lower()      # LOW confidence (lowercase)
    if inner == "N" and outer == "N":
        return "N"                # Both missing
    return "n"                    # Disagreement (unresolved)
```

### Worked example: full pipeline for one position

Let us trace position chr1:50,000 through the entire algorithm:

```text
  Step 1: Collect projected bases at chr1:50,000

    Bonobo (inner):    T
    Chimp (inner):     T
    Gorilla (inner):   T
    Macaque (outer):   T

  Step 2: Inner majority vote

    Tally: T=3
    min_inner_freq=1, so T qualifies
    Inner consensus: T

  Step 3: Outer majority vote

    Tally: T=1
    min_outer_freq=1, so T qualifies
    Outer consensus: T

  Step 4: Compare

    Inner=T, Outer=T → they agree
    Output: "T" (uppercase = high confidence)
```

Now consider a more interesting position, chr1:75,000:

```text
  Step 1: Collect projected bases at chr1:75,000

    Bonobo (inner):    A
    Chimp (inner):     A
    Gorilla (inner):   G
    Macaque (outer):   G

  Step 2: Inner majority vote

    Tally: A=2, G=1
    Inner consensus: A (majority)

  Step 3: Outer majority vote

    Tally: G=1
    Outer consensus: G

  Step 4: Compare

    Inner=A, Outer=G → DISAGREEMENT
    Output: "n" (unresolved)
```

This pattern (inner says A, outer says G) is consistent with incomplete lineage sorting in the human-gorilla ancestor, where gorilla and macaque inherited one allele and bonobo+chimp inherited another.

---

## Confidence encoding

The output FASTA encodes confidence via letter case. This is compatible with the Ensembl EPO ancestral sequence convention:

| Character | Confidence | Condition | Typical fraction (human) |
|-----------|-----------|-----------|--------------------------|
| `ACGT` | **High** | Inner and outer outgroups agree | ~75% |
| `acgt` | **Low** | Only one outgroup tier has data | ~15% |
| `n` | **Unresolved** | Inner and outer disagree | ~0.1% |
| `N` | **Missing** | Both tiers lack data | ~10% |

This encoding means that a single FASTA file carries both the ancestral call and its reliability. No sidecar files or separate quality tracks are needed.

---

## Biological rationale

### Why two tiers instead of one big vote?

Consider what happens with a single flat vote across all outgroups:

```text
  Single-tier vote (all 4 species equally weighted):

    Bonobo:   A     ┐
    Chimp:    A     ├─ These 3 are NOT independent observations.
    Gorilla:  A     ┘  Due to ILS, they may all inherit the same
    Macaque:  G        wrong allele from the ancestral population.

    Flat vote: A wins 3-to-1.
    But what if A is the derived allele, shared through ILS?
```

The two-tier approach treats the inner outgroups as a *single correlated observation* and demands independent confirmation from a distant outgroup:

```text
  Two-tier vote:
    Inner consensus: A (3/3)
    Outer consensus: G (1/1)
    → DISAGREEMENT → flagged as "n"

  This catches the ILS case!
```

### Incomplete lineage sorting (ILS)

When ancestral populations were large and speciation events were close together, gene trees can disagree with the species tree:

```text
  Species tree:                  Gene tree at this locus:

       ┌── Bonobo                     ┌── Bonobo
    ┌──┤                           ┌──┤
    │  └── Chimp                   │  └── Gorilla     ← ILS!
  ──┤                            ──┤
    │  ┌── Gorilla                 │  ┌── Chimp
    └──┤                           └──┤
       └── Macaque                    └── Macaque
```

For the human-chimp-gorilla clade, ILS affects ~1-3% of the genome. The two-tier approach detects most of these cases because the outer outgroup (macaque) provides a phylogenetically independent anchor.

### When the algorithm can still fail

No method is perfect. The algorithm can produce incorrect calls in rare cases:

| Failure mode | Frequency | Mechanism |
|-------------|-----------|-----------|
| Convergent substitution | <0.01% | Same mutation on both inner and outer branches |
| Correlated ILS | <0.1% | All inner species AND the outer share the derived allele through deep coalescence |
| Systematic alignment error | Variable | Paralog confusion in segmental duplications |

These contribute to the ~0.1% disagreement rate observed when comparing ancify's human calls to the Ensembl EPO 13-primate reference.

---

## The `min_inner_freq` parameter in depth

This parameter lets you trade coverage for accuracy:

```text
  3 inner species, varying min_inner_freq:

  ┌───────────────────┬───────────┬──────────┬──────────────────────┐
  │ min_inner_freq    │ Rule      │ Coverage │ Accuracy             │
  ├───────────────────┼───────────┼──────────┼──────────────────────┤
  │ 1                 │ Any 1     │ Highest  │ Lowest (but still    │
  │                   │ suffices  │          │ >99% for primates)   │
  │ 2                 │ 2 of 3    │ Moderate │ Good                 │
  │                   │ agree     │          │                      │
  │ 3                 │ All 3     │ Lowest   │ Highest              │
  │                   │ agree     │          │                      │
  └───────────────────┴───────────┴──────────┴──────────────────────┘
```

**Example:** With `min_inner_freq=2` and bases `[A, N, A]`:
- Tally: A=2, N excluded
- Threshold met (2 >= 2): consensus is A

With `min_inner_freq=2` and bases `[A, N, G]`:
- Tally: A=1, G=1
- Neither reaches threshold: consensus is N (missing)

The same logic applies to `min_outer_freq` for the outer tier. Since there is typically only one outer outgroup species, `min_outer_freq` is usually left at 1.

---

## Alignment quality and its effects

The pipeline's accuracy depends fundamentally on the quality of the input alignments:

```text
  Good alignment region        Poor alignment region
  (unique sequence)            (segmental duplication)

  Focal:    ACGTACGT            Focal:    ACGTACGT
  Outgroup: ACGTACGT            Outgroup: TGCATGCA  ← paralog!
            ^^^^^^^^                      ^^^^^^^^
            Correct projected             Incorrect projected
            bases                         bases → bad ancestral call
```

Net alignments mitigate paralog confusion by keeping only the best one-to-one chain at each position, but regions with recent segmental duplications, transposable elements, or structural variants can still produce artifacts. These typically manifest as `N` (no alignment) or disagreements (`n`).

---

## Alternative method: Fitch parsimony

Since version 1.2.0 ancify supports a second inference method: **Fitch parsimony** (Fitch, 1971). Instead of splitting outgroups into two hand-picked tiers and majority-voting within each, this method uses the actual **phylogenetic tree** of the outgroup species to reconstruct the most parsimonious ancestral allele at the root.

### When to prefer parsimony

- You have **three or more outgroups** at varying phylogenetic distances.
- You want the tree topology to determine how species are weighted, rather than manually choosing inner/outer groups.
- Resolving cases that the two-tier method marks as "unresolved" (`n`) matters to you. In many of those cases, parsimony can find a unique most-parsimonious solution by leveraging the tree structure.

### The Fitch algorithm

The algorithm runs on a rooted phylogenetic tree whose leaves are the outgroup species. At every genomic position it performs two passes:

**Pass 1: Bottom-up (post-order traversal)**

Starting at the leaves, work toward the root. Each leaf gets a set containing its observed allele (or `{A,C,G,T}` if missing). At each internal node, compute:

- **Intersection** of children's sets, if non-empty
- **Union** of children's sets, otherwise

The intersection case means all children are compatible; no mutation is needed on the branch. The union case means at least one mutation is required.

```text
  Example: (((bonobo,chimp),gorilla),macaque)
  Observed: bonobo=G, chimp=G, gorilla=A, macaque=A

  Bottom-up:
    bonobo  → {G}       chimp → {G}
    (bonobo,chimp)      → {G} ∩ {G} = {G}
    gorilla → {A}
    ((bc),gorilla)      → {G} ∩ {A} = ∅ → {G} ∪ {A} = {A,G}
    macaque → {A}
    root                → {A,G} ∩ {A} = {A}
```

**Pass 2: Top-down (pre-order traversal)**

Starting at the root, assign a concrete allele from the node's set. Prefer the parent's allele if it is in the child's set (minimizes mutations). Ties are broken alphabetically for determinism.

```text
  Continuing the example:
    root  → {A}    → picks A
    (bc,gorilla) → {A,G} → parent is A, A ∈ set → picks A
    gorilla      → {A}   → A
    (bc)         → {G}   → parent is A, A ∉ set → picks G
    bonobo       → {G}   → G
    chimp        → {G}   → G

  Result: A is ancestral at the root (1 mutation: A→G on the bonobo-chimp branch)
```

### Handling missing data

When a leaf has no alignment data (`N`), it is assigned the full set `{A,C,G,T}`. This makes it compatible with any allele, so missing data never forces a mutation. If **all** leaves are missing, the root set remains `{A,C,G,T}` and the position is reported as `N` (missing).

### Confidence encoding (parsimony)

The parsimony method maps to the same case-based confidence encoding as the voting method:

| Character | Confidence | Condition |
|-----------|-----------|-----------|
| `ACGT` | **High** | Root Fitch set has exactly one allele (unambiguous) |
| `acgt` | **Low** | Root Fitch set has multiple alleles (ambiguous; one is chosen) |
| `N` | **Missing** | All outgroup leaves are `N` |

Note that there is no `n` (unresolved) category for parsimony — the algorithm always produces a call unless all data is missing.

### Comparison: voting vs. parsimony

```text
  Position with: bonobo=G, chimp=G, gorilla=A, macaque=A

  Two-tier voting:
    Inner consensus (bonobo, chimp, gorilla) → majority = G
    Outer consensus (macaque)                → A
    Inner ≠ Outer → "n" (UNRESOLVED)

  Fitch parsimony on (((bonobo,chimp),gorilla),macaque):
    Root = A (one mutation needed: A→G on the bonobo-chimp branch)
    → "A" (HIGH CONFIDENCE)
```

The tree structure lets parsimony recognize that the G allele is a shared derived change in the bonobo-chimp clade, while A is the ancestral state supported by both gorilla and macaque.

### Configuration

To use parsimony, add two fields to your YAML config:

```yaml
method: parsimony
tree: "(((bonobo,chimp),gorilla),macaque)"
```

The `tree` field accepts either an inline Newick string or a path to a `.nwk` file. Leaf names must match outgroup `name` fields in your config. See {doc}`configuration` for the complete reference.

---

## Summary

### Two-tier voting (default)

```text
  For each position in the focal genome:

  1. Gather projected bases from all outgroup species
  2. Compute inner consensus (majority vote among close relatives)
  3. Compute outer consensus (majority vote among distant relatives)
  4. Compare:
       agree  → UPPERCASE (high confidence)
       one N  → lowercase (low confidence)
       differ → n (unresolved)
       both N → N (missing)
  5. Write to output FASTA
```

### Fitch parsimony

```text
  For each position in the focal genome:

  1. Gather projected bases from all outgroup species
  2. Run Fitch bottom-up pass on the phylogenetic tree
  3. Run Fitch top-down pass to assign the root allele
  4. Encode confidence:
       root set size = 1 → UPPERCASE (high confidence)
       root set size > 1 → lowercase (low confidence)
       all missing       → N
  5. Write to output FASTA
```

Both methods are simple, fast (Phase 2 takes minutes for a full genome), and produce reliable ancestral calls with built-in confidence assessment. Parsimony is more principled for complex outgroup configurations; voting is simpler and well-tested for the common two-tier setup.
