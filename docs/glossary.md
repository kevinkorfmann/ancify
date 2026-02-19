# Glossary

A reference for terms used throughout the ancify documentation.

---

Ancestral allele
: The nucleotide present at a genomic position in the common ancestor of the species being studied. Inferred by comparison with outgroup species.

AXT format
: A text-based pairwise alignment format used by UCSC. Each alignment block has a header line followed by the target sequence, the query sequence, and a blank line. ancify reads **net AXT** files, which represent best-in-genome one-to-one alignments.

Confidence encoding
: ancify's convention for encoding both the ancestral call and its reliability in a single character. Uppercase (`ACGT`) = high confidence, lowercase (`acgt`) = low confidence, `n` = unresolved, `N` = missing.

Convergent substitution
: An identical mutation arising independently on two separate lineages. When this happens on both the inner and outer outgroup branches, it can cause a false high-confidence call. Rare in practice (<0.01% of sites).

Coordinate projection
: The process of mapping an outgroup species' aligned bases onto the focal species' genomic coordinates (Phase 1). The result is a FASTA file with the same length as the focal chromosome, where each position contains the outgroup's base (or `N` if there is no alignment).

Derived allele
: A nucleotide that arose by mutation since the common ancestor. The complement of the ancestral allele.

Focal species
: The species whose genome you are polarizing. All coordinates and chromosome names refer to this species' reference assembly.

Folded SFS
: The site frequency spectrum where minor and major allele frequencies are combined (symmetric around 0.5). Does not require polarization, but carries less information than the unfolded SFS.

ILS (Incomplete Lineage Sorting)
: A phenomenon where gene trees differ from the species tree due to ancestral polymorphism persisting across speciation events. More common when ancestral population sizes are large and inter-speciation times are short. ancify detects potential ILS by flagging inner-outer disagreements as `n`.

Inner outgroup
: One or more species closely related to the focal species (e.g. bonobo, chimp, gorilla for a human focal species). Their bases are combined by majority vote to form the inner consensus.

Majority vote
: The algorithm used to combine bases from multiple species in the same outgroup tier. The most frequent valid nucleotide (A, C, G, T) wins. Ties are broken alphabetically. Bases below `min_inner_freq` or `min_outer_freq` are treated as missing.

`min_inner_freq` / `min_outer_freq`
: Configuration parameters controlling the minimum allele count required for a majority vote to produce a consensus. Higher values require more species to agree, increasing stringency but reducing coverage.

Net alignment
: A processed alignment where overlapping chains have been resolved to produce one-to-one best matches. "Net" refers to the netting algorithm that removes lower-scoring overlaps. Net AXT files from UCSC are the standard input for ancify.

Outer outgroup
: One or more species distantly related to the focal species (e.g. macaque for a human focal species). Provides an independent evolutionary check on the inner consensus.

Parsimony
: The principle that the simplest explanation is preferred. In ancestral allele inference, if multiple outgroup species share the same allele, parsimony infers that this was the ancestral state.

Phase 1 / Phase 2 / Phase 3
: The three stages of the ancify pipeline. **Phase 1** (project): coordinate projection of outgroup alignments. **Phase 2** (call): ancestral state inference via two-tier voting. **Phase 3** (evaluate): optional comparison against external references.

Polarization
: The process of determining which allele at a variant site is ancestral and which is derived. Required for computing the unfolded site frequency spectrum.

Site Frequency Spectrum (SFS)
: A summary statistic counting the number of variant sites at each allele frequency in a sample. The **unfolded** SFS uses derived allele frequencies and requires polarization. Used extensively in demographic inference and tests of selection.

Two-tier voting
: ancify's approach of splitting outgroup species into two independent groups (inner and outer), computing a majority vote within each, and comparing the results. Agreement yields high confidence; disagreement yields an unresolved flag.

Unfolded SFS
: The site frequency spectrum expressed in terms of **derived** allele frequencies. Requires knowing which allele is ancestral (polarization). More informative than the folded SFS for demographic inference and selection analysis.

UCSC Genome Browser
: A web-based genomic data platform ([genome.ucsc.edu](https://genome.ucsc.edu/)) that hosts genome assemblies, annotations, and pairwise alignments for hundreds of species. The primary source of net AXT alignment files used by ancify.
