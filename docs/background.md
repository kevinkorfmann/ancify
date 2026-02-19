# Population Genetics Background

*This page teaches you the biology behind ancestral allele polarization. If you already know what the site frequency spectrum is and why it needs to be unfolded, feel free to skip to the {doc}`quickstart`.*

---

## The big picture

Every individual in a population carries DNA that differs — by single nucleotide changes, insertions, deletions, and rearrangements — from every other individual. Population geneticists study these **variants** to learn about evolution, natural selection, demographic history, and disease.

But there is a question that comes up immediately: when we see a variant, **which allele is the original one (ancestral) and which is the new mutation (derived)?**

This matters because the frequency of derived alleles in a population is profoundly informative. A derived allele at 1% frequency tells a very different story than one at 99% — but if we cannot distinguish ancestral from derived, we cannot tell these apart.

## Alleles, mutations, and time

Consider a position in the human genome where some people carry **A** and others carry **G**:

```text
Population sample at position chr1:1,000,000

  Individual 1:   A
  Individual 2:   A
  Individual 3:   G    ← this one is different
  Individual 4:   A
  Individual 5:   A
```

One of two things happened:

1. The **ancestral** allele was **A**, and a mutation created **G** (G is *derived*).
2. The **ancestral** allele was **G**, and a mutation created **A** (A is *derived*).

These two scenarios lead to very different conclusions about whether this variant is under selection, how old it is, and what it tells us about population history.

## The site frequency spectrum (SFS)

The **site frequency spectrum** counts how many variants exist at each frequency in the population. It is one of the most important summaries in population genetics.

```text
                        Folded SFS
       (cannot distinguish rare derived from rare ancestral)

 Count │
  of   │ ████
 sites │ ████ ███
       │ ████ ███ ██
       │ ████ ███ ██  █    █
       └──────────────────────────
         0.1  0.2  0.3  0.4  0.5
              Minor allele frequency


                       Unfolded SFS
           (polarized: we know which allele is derived)

 Count │
  of   │ ████
 sites │ ████ ███
       │ ████ ███ ██
       │ ████ ███ ██  █    █         █
       └──────────────────────────────────
         0.1  0.2  0.3  0.4  0.5 ... 0.9
              Derived allele frequency
```

The **folded** SFS (top) is symmetric around 0.5 — you can always compute it. But the **unfolded** SFS (bottom) requires knowing which allele is ancestral. It carries much more information:

- **Neutral theory predicts** that the unfolded SFS follows a 1/f distribution (many rare derived alleles, few common ones).
- **Positive selection** creates an excess of high-frequency derived alleles.
- **Purifying selection** removes derived alleles, depleting the intermediate frequencies.
- **Demographic history** (bottlenecks, expansions, migration) leaves distinctive signatures in the SFS shape.

Tools like **dadi**, **moments**, **fastsimcoal2**, and **mushi** all require the *unfolded* SFS — which means they need polarization.

## How do we determine the ancestral allele?

We cannot go back in time, but we can use **outgroup species** — organisms that share a common ancestor with our focal species but diverged long ago.

### The parsimony argument

```text
                     ┌─── Human (A/G polymorphic)
           ┌────────┤
           │        └─── Chimp (A)
     ──────┤
           │        ┌─── Gorilla (A)
           └────────┤
                    └─── Macaque (A)

  If chimp, gorilla, and macaque all carry A,
  parsimony says A is the ancestral allele.
  The G in some humans is derived.
```

The logic: if multiple outgroup species all carry the same allele, the simplest explanation is that this was the ancestral state, and the different allele in the focal species arose by mutation.

### Why multiple outgroups?

Using a single outgroup is risky:

| Problem | What happens | How multiple outgroups help |
|---------|-------------|---------------------------|
| **Back mutation** | The outgroup mutated away from the true ancestral allele | Other outgroups still carry the correct allele |
| **Incomplete lineage sorting (ILS)** | The gene tree disagrees with the species tree | Majority vote across several species reduces error |
| **Alignment error** | The outgroup sequence at this position is wrong | Redundancy from other species compensates |

### The two-tier approach

ancify goes a step further by splitting outgroups into two tiers:

```text
  INNER outgroups (close relatives)     OUTER outgroup (distant relative)
  ──────────────────────────────────    ────────────────────────────────
  Bonobo:  A  ─┐                        Macaque: A
  Chimp:   A  ─┼─▶ majority vote: A           │
  Gorilla: A  ─┘                               │
                                                │
            ┌───────────────────────────────────┘
            ▼
     Inner says A, Outer says A → Agreement!
     → High-confidence ancestral call: A (uppercase)
```

**Why two tiers?** The inner outgroups are close enough that their alignments cover most of the genome, but close enough that ILS can affect them as a group. The outer outgroup is far enough away that ILS between the tiers is negligible — so when inner and outer *agree*, you can be very confident.

When they *disagree*, something interesting is happening (ILS, convergent mutation, alignment artifact), and the position is flagged as unresolved rather than guessed at.

## What about non-model organisms?

The same logic works for *any* species with available outgroup alignments:

| Focal species | Inner outgroups | Outer outgroup |
|--------------|----------------|----------------|
| Human | Bonobo, Chimp, Gorilla (~6-9 Mya) | Macaque (~25 Mya) |
| Mouse | Rat (~12 Mya) | Rabbit (~90 Mya) |
| *D. melanogaster* | *D. simulans*, *D. sechellia* (~2.5 Mya) | *D. yakuba* (~6 Mya) |
| *Brassica rapa* | *B. oleracea* (same genus) | *Arabidopsis thaliana* (~20 Mya) |
| Zebrafish | Medaka | Fugu |

The UCSC Genome Browser provides pairwise alignments for hundreds of species pairs. If your focal species has a UCSC assembly, you likely already have the input data you need.

## The confidence encoding

ancify encodes both the ancestral call *and* its confidence in a single character:

```text
  Position 1:  Inner=A, Outer=A  →  "A"  (both agree: HIGH confidence)
  Position 2:  Inner=A, Outer=N  →  "a"  (outer missing: LOW confidence)
  Position 3:  Inner=N, Outer=T  →  "t"  (inner missing: LOW confidence)
  Position 4:  Inner=A, Outer=T  →  "n"  (disagree: UNRESOLVED)
  Position 5:  Inner=N, Outer=N  →  "N"  (both missing: NO DATA)
```

This means you can choose your stringency downstream:

- **Conservative analysis**: use only uppercase positions (inner+outer agree).
- **Maximum coverage**: include lowercase positions (one-tier data).
- **Exclude trouble spots**: skip `n` positions (disagreements flag potential ILS or alignment problems).

## Where does ancify fit in a typical workflow?

```text
  1. Obtain genome assembly          (e.g. hg38 from UCSC)
  2. Download outgroup alignments    (net AXT files from UCSC)
  3. Run ancify                      ← you are here
  4. Get ancestral FASTA files       (one per chromosome)
  5. Polarize your VCF variants      (annotate REF/ALT as ancestral/derived)
  6. Compute unfolded SFS            (input for demographic inference)
  7. Run downstream analyses         (dadi, moments, selection scans, etc.)
```

ancify handles step 3, transforming raw alignment data into a per-position ancestral reference that plugs directly into your existing variant-calling and analysis pipelines.

## Ready to start?

Head to the {doc}`quickstart` to go from zero to ancestral alleles in five minutes, or dive into the {doc}`algorithm` for the full technical details.
