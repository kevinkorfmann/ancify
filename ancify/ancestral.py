"""Phase 2: Infer ancestral alleles from projected outgroup sequences.

Uses a two-tier outgroup voting scheme:

* **Inner outgroup** -- closely related species; the most frequent
  nucleotide among them forms the inner consensus.
* **Outer outgroup** -- more distantly related species; serves as an
  independent confirmation.

Confidence is encoded via letter case in the output FASTA:

========  ===========  ==========================================
Char      Confidence   Condition
========  ===========  ==========================================
``ACGT``  High         Inner and outer outgroups agree
``acgt``  Low          Only one tier has data
``n``     Unresolved   Inner and outer disagree
``N``     Missing      Both tiers lack data
========  ===========  ==========================================
"""

import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from .utils import read_fasta, write_fasta, majority_vote, VALID_ALLELES

logger = logging.getLogger(__name__)


def call_ancestral_base(inner_bases, outer_bases,
                        min_inner_freq=1, min_outer_freq=1):
    """Infer the ancestral allele at a single position.

    Parameters
    ----------
    inner_bases : list of str
        Nucleotides from the inner (closely related) outgroup species.
    outer_bases : list of str
        Nucleotides from the outer (distantly related) outgroup species.
    min_inner_freq, min_outer_freq : int
        Minimum allele count to accept a majority-vote consensus.

    Returns
    -------
    str
        Single character with case-encoded confidence (see module docstring).
    """
    inner = majority_vote(inner_bases, min_inner_freq)
    outer = majority_vote(outer_bases, min_outer_freq)

    if inner != "N" and inner == outer:
        return inner
    if inner == "N" and outer != "N":
        return outer.lower()
    if inner != "N" and outer == "N":
        return inner.lower()
    if inner == "N" and outer == "N":
        return "N"
    return "n"


def _call_chromosome(args):
    """Worker: call ancestral states for one chromosome."""
    chrom, inner_paths, outer_paths, out_path, min_inner, min_outer = args

    inner_seqs = [read_fasta(p)[1] for p in inner_paths]
    outer_seqs = [read_fasta(p)[1] for p in outer_paths]

    length = len(inner_seqs[0])
    for s in inner_seqs + outer_seqs:
        if len(s) != length:
            raise ValueError(
                f"Length mismatch on {chrom}: expected {length}, got {len(s)}"
            )

    anc = []
    for i in range(length):
        ib = [s[i].upper() for s in inner_seqs]
        ob = [s[i].upper() for s in outer_seqs]
        anc.append(call_ancestral_base(ib, ob, min_inner, min_outer))

    write_fasta(out_path, f">{chrom}", "".join(anc))
    return chrom


def run_ancestral_calling(config):
    """Execute Phase 2: call ancestral alleles for every chromosome.

    Reads projected FASTA files from ``<work_dir>/projected/<species>/``
    and writes ancestral FASTA files to ``<output_dir>/<chrom>.fa``.
    """
    chromosomes = config.resolve_chromosomes()
    work = Path(config.work_dir)
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    for chrom in chromosomes:
        inner_paths = [
            str(work / "projected" / og.name / f"{chrom}.fa")
            for og in config.outgroups_inner
        ]
        outer_paths = [
            str(work / "projected" / og.name / f"{chrom}.fa")
            for og in config.outgroups_outer
        ]
        out_path = str(out_dir / f"{chrom}.fa")
        tasks.append((
            chrom, inner_paths, outer_paths, out_path,
            config.min_inner_freq, config.min_outer_freq,
        ))

    logger.info("Phase 2: calling ancestral states for %d chromosomes", len(tasks))

    with ProcessPoolExecutor(max_workers=config.num_cpus) as pool:
        futures = {pool.submit(_call_chromosome, t): t for t in tasks}
        for future in as_completed(futures):
            chrom = future.result()
            logger.info("  Completed %s", chrom)

    logger.info("Phase 2 complete.")
