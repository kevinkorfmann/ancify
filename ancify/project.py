"""Phase 1: Project outgroup alignments onto focal-species coordinates.

Reads pairwise net-AXT alignment files and produces one FASTA file per
chromosome, where each position corresponds to a position in the focal
species reference.  Unaligned positions are filled with ``N``.
"""

import gzip
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from .utils import read_chromosome_lengths

logger = logging.getLogger(__name__)


def project_alignment(axt_path, chrom_length, target_chrom):
    """Project one chromosome from an AXT alignment onto focal coordinates.

    Parameters
    ----------
    axt_path : str or Path
        Path to the (optionally gzip-compressed) net-AXT alignment file.
    chrom_length : int
        Length of *target_chrom* in the focal species reference.
    target_chrom : str
        Chromosome name to extract (must match the target field in the AXT).

    Returns
    -------
    str
        Outgroup sequence in focal-species coordinates (length == *chrom_length*).
        Positions not covered by any alignment block contain ``N``.
    """
    seq = bytearray(b"N" * chrom_length)

    opener = gzip.open if str(axt_path).endswith(".gz") else open
    linecount = 0
    hc, hs = None, 0
    hseq, cseq = None, None

    with opener(axt_path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue

            phase = linecount % 4
            if phase == 0:
                parts = line.split()
                hc = parts[1]
                hs = int(parts[2])
            elif phase == 1:
                hseq = line
            elif phase == 2:
                cseq = line
            elif phase == 3:
                if hc == target_chrom and hseq and cseq:
                    pos = hs
                    for i in range(len(hseq)):
                        if hseq[i].upper() in "ATCGN":
                            idx = pos - 1
                            if 0 <= idx < chrom_length:
                                seq[idx] = ord(cseq[i])
                            pos += 1
                hc, hs = None, 0
                hseq, cseq = None, None
            linecount += 1

    return seq.decode("ascii")


def _project_one(args):
    """Worker: project a single chromosome for a single species."""
    axt_path, chrom_length, chrom, out_path = args
    logger.info("Projecting %s from %s", chrom, Path(axt_path).name)
    sequence = project_alignment(axt_path, chrom_length, chrom)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f">{chrom}\n{sequence}\n")
    return chrom


def run_projection(config):
    """Execute Phase 1 for all outgroup species and chromosomes.

    Creates ``<work_dir>/projected/<species_name>/<chrom>.fa`` for every
    combination of outgroup species and chromosome.
    """
    chrom_lengths = read_chromosome_lengths(config.chromosome_lengths)
    chromosomes = config.resolve_chromosomes()
    work = Path(config.work_dir)

    all_outgroups = config.outgroups_inner + config.outgroups_outer
    tasks = []

    for og in all_outgroups:
        species_dir = work / "projected" / og.name
        species_dir.mkdir(parents=True, exist_ok=True)
        for chrom in chromosomes:
            if chrom not in chrom_lengths:
                logger.warning("Chromosome %s not in lengths file, skipping", chrom)
                continue
            out_path = species_dir / f"{chrom}.fa"
            tasks.append((og.alignment, chrom_lengths[chrom], chrom, str(out_path)))

    n_species = len(all_outgroups)
    n_chroms = len(chromosomes)
    logger.info(
        "Phase 1: projecting %d tasks (%d species x %d chromosomes)",
        len(tasks), n_species, n_chroms,
    )

    with ProcessPoolExecutor(max_workers=config.num_cpus) as pool:
        futures = {pool.submit(_project_one, t): t for t in tasks}
        for future in as_completed(futures):
            chrom = future.result()
            logger.info("  Completed %s", chrom)

    logger.info("Phase 1 complete.")
