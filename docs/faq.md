# FAQ & Troubleshooting

Common questions and solutions for issues you might encounter.

---

## General Questions

### What species can I use ancify with?

Any species for which you have:
1. A reference genome assembly
2. Pairwise net AXT alignments to at least two outgroup species (one inner, one outer)

The UCSC Genome Browser provides alignments for hundreds of species pairs. If your focal species has a UCSC assembly, you almost certainly have the data you need.

### How accurate is ancify?

For the human genome (hg38) using the BCGM method (bonobo + chimp + gorilla + macaque):

- **>99.9% agreement** with the Ensembl EPO 13-primate ancestral reference at high-confidence sites
- **~0.08% disagreement rate** on autosomes
- **>99.5%** of ancestral calls match either the REF or ALT allele at known variant sites

Accuracy depends on the divergence times and alignment quality for your specific species.

### How long does it take?

| Phase | CPU (scalar) | CPU (vectorized) | GPU (NVIDIA A100) |
|-------|-------------|-----------------|-------------------|
| Phase 1 (projection) | 2-8 hours | ~5 minutes | ~2 minutes |
| Phase 2 (calling) | 5-15 minutes | ~10 minutes | ~10 seconds |
| Phase 3 (evaluation) | 5-15 minutes | 5-15 minutes | 5-15 minutes |

With the GPU backend enabled, the full human genome completes in ~2 minutes. Even without a GPU, the vectorized NumPy backend is much faster than the original scalar path. See {doc}`performance` for details.

### How much memory do I need?

Phase 2 is the memory-intensive step:

```text
memory ≈ num_cpus × (num_outgroups) × (avg_chrom_length)
```

For human chr1 (~249 Mb) with 4 outgroups and 24 CPUs: ~24 GB peak. Reduce `num_cpus` if needed.

### Can I run phases independently?

Yes. Each phase has its own subcommand:

```bash
ancify project  -c config.yaml   # Phase 1 only
ancify call     -c config.yaml   # Phase 2 only (requires Phase 1 output)
ancify evaluate -c config.yaml   # Phase 3 only (requires Phase 2 output)
```

This is useful when you want to re-run Phase 2 with different `min_inner_freq` settings without repeating the expensive Phase 1.

---

## Data Questions

### Where do I get net AXT alignment files?

From the UCSC Genome Browser downloads:

```text
https://hgdownload.soe.ucsc.edu/goldenPath/<focal_assembly>/vs<Outgroup>/
```

For example, human vs. chimp:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro6/hg38.panTro6.net.axt.gz
```

Browse all available alignments at [UCSC Downloads](https://hgdownload.soe.ucsc.edu/downloads.html).

### What if UCSC does not have alignments for my species?

You can generate your own net AXT files using the UCSC Kent tools:

1. Run `lastz` or `minimap2` to produce pairwise alignments
2. Convert to chain format with `axtChain`
3. Net with `chainNet`
4. Convert back to AXT with `netToAxt`

This is more involved but allows ancify to work with any pair of genome assemblies.

### Can I use MAF files instead of AXT?

Not directly. ancify reads AXT format. You can convert MAF to AXT using UCSC tools:

```bash
mafToAxt input.maf focal_assembly outgroup_assembly output.axt
```

### What chromosome naming convention should I use?

ancify does not enforce any naming convention. Use whatever names appear in your chromosome lengths file. Examples:

- `chr1`, `chr2`, ..., `chrX` (UCSC style)
- `1`, `2`, ..., `X` (Ensembl style)
- `2L`, `2R`, `3L`, `3R` (Drosophila)
- `A01`, `A02` (Brassica)

The names just need to be consistent between your lengths file, your alignment files, and your config.

---

## Performance & Resource Issues

### `MemoryError` during Phase 2

Each chromosome loads all projected sequences simultaneously. Solutions:

1. **Reduce `num_cpus`** — this limits how many chromosomes are loaded in parallel.
2. **Process a subset of chromosomes** — use the `chromosomes` config field to run in batches.
3. **Use a machine with more RAM** — Phase 2 for human chr1 with 4 outgroups needs ~1 GB per concurrent chromosome.

### Phase 1 is very slow

Phase 1 streams through large compressed files. This is I/O-bound, not CPU-bound:

- **Install `isal`** — `pip install 'ancify[fast]'` provides 2–5× faster gzip decompression using Intel ISA-L
- **Enable the vectorized backend** — set `backend: auto` or `backend: cpu` in your config. The vectorized NumPy scatter is 20–50× faster than the original character loop
- Use an **SSD** for the alignment files
- Consider splitting by species: run projection for each outgroup on a separate machine/job, since AXT files are independent
- Phase 1 only needs to be run once; subsequent Phase 2 runs reuse the projected files

### Do I need a GPU?

No. The vectorized NumPy backend (selected automatically when PyTorch is not installed, or with `backend: cpu`) is already much faster than the original scalar path. GPU acceleration is an additional speedup on top of that, primarily beneficial for large genomes (> 1 Gb). See {doc}`performance` for benchmarks.

### How do I enable GPU acceleration?

1. Install PyTorch with CUDA support: `pip install torch --index-url https://download.pytorch.org/whl/cu128`
2. Set `backend: auto` (or `backend: gpu`) in your config
3. Optionally set `gpu_devices: [0, 1, 2]` to specify which GPUs to use

ancify will detect PyTorch + CUDA automatically and use the GPU backend. Check with `ancify run -c config.yaml -v` — the first log line reports the active backend.

### My GPU runs out of memory

Each chromosome needs roughly 6 GB of GPU memory for 4 outgroups and a 248M-position chromosome. Solutions:

1. **Process fewer chromosomes in parallel** — reduce `num_cpus`
2. **Restrict to fewer GPUs** — set `gpu_devices` to a list with fewer entries so memory is not split across too many concurrent tasks
3. **Fall back to CPU** — set `backend: cpu` for the vectorized NumPy path, which uses system RAM instead

### Does the GPU backend change the output?

No. All three backends (scalar, vectorized CPU, GPU) produce **bit-identical** output. The tie-breaking rule, frequency thresholds, and confidence encoding are preserved exactly. You can verify this by running the same config with different `backend` settings and comparing the output FASTAs.

### Can I resume a failed run?

ancify does not have built-in checkpointing, but:

- **Phase 1** writes one file per (species, chromosome). If some files already exist, you can manually skip those species by temporarily removing them from the config and re-adding after.
- **Phase 2** is fast enough to re-run from scratch.
- **Phase 3** is also fast and can be re-run independently.

---

## Output Questions

### What does each character mean in the output FASTA?

| Character | Confidence | Condition |
|-----------|-----------|-----------|
| `A` `C` `G` `T` | **High** | Inner and outer outgroups agree |
| `a` `c` `g` `t` | **Low** | Only one outgroup tier has data |
| `n` | **Unresolved** | Inner and outer disagree |
| `N` | **Missing** | Both tiers lack data |

### My output has a lot of `N`s. Is something wrong?

Not necessarily. High `N` rates are expected in:

- **Repetitive regions** (centromeres, telomeres, satellite DNA)
- **Lineage-specific insertions** (transposable elements inserted in the focal species after divergence from outgroups)
- **Assembly gaps** in the focal genome
- **chrY** (very repetitive, poor alignment coverage)
- **Species with distant outgroups** (less alignable sequence)

Compare your coverage rates against the expected values for your species pair's divergence time.

### Can I use the output to polarize a VCF?

Yes — this is the primary use case. See the {doc}`quickstart` for code examples and the {doc}`tutorials` for a complete walkthrough.

### Is the confidence encoding compatible with Ensembl EPO?

Yes. The uppercase/lowercase/N convention matches the Ensembl EPO ancestral sequence format. Tools that parse Ensembl ancestral sequences should work with ancify output without modification.

---

## Scientific Questions

### How does this compare to the Ensembl EPO method?

The Ensembl EPO (Enredo-Pecan-Ortheus) method uses:
- Multiple sequence alignment of 13+ primate genomes
- A probabilistic model (Ortheus) that reconstructs ancestral sequences on a phylogenetic tree
- An order of magnitude more species

ancify uses:
- Pairwise alignments only (simpler input, widely available)
- Deterministic majority voting (transparent, reproducible)
- Typically 2-5 outgroup species

**Trade-offs:**

| Aspect | Ensembl EPO | ancify |
|--------|------------|--------|
| Input complexity | High (MSA of 13+ genomes) | Low (2-5 pairwise AXTs) |
| Species coverage | Vertebrates only | Any species with AXT alignments |
| Accuracy | Slightly higher | Very close (>99.9% agreement) |
| Transparency | Probabilistic model | Simple, auditable rules |
| Customizability | Fixed | Full control via config |

For most applications, the practical difference in accuracy is negligible.

### When should I NOT use parsimony-based polarization?

Parsimony can be unreliable when:

1. **Divergence times are very large** — back mutations accumulate, eroding the ancestral signal.
2. **The focal species is in a rapid radiation** — ILS is pervasive and multiple outgroups may all carry derived alleles.
3. **Only one outgroup is available and it is distant** — error rates increase without redundancy.

In these cases, consider probabilistic methods (like Ortheus or est-sfs) that explicitly model substitution rates along branches.

### What `min_inner_freq` should I use?

It depends on your tolerance for errors vs. missing data:

- **min_inner_freq=1** (default): Maximize coverage. Best for demographic inference where the unfolded SFS shape matters more than individual-site accuracy.
- **min_inner_freq=2** (with 3 inner species): A good balance. Recommended if you are doing site-by-site analysis (e.g. McDonald-Kreitman tests).
- **min_inner_freq=N** (unanimity): Maximum stringency. Use when even rare errors would be problematic.

---

## Still stuck?

Open an issue on [GitHub](https://github.com/kevinkorfmann/ancify/issues) with:
1. Your config file (with paths anonymized if needed)
2. The full error message
3. The output of `ancify --help` (to confirm installation)
4. Your Python version (`python --version`)
