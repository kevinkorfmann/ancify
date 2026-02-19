"""I/O utilities and common helpers for ancify."""

import gzip
from pathlib import Path

VALID_ALLELES = ("A", "C", "G", "T")


def read_fasta(path):
    """Read a single-record FASTA file and return (header, sequence).

    Handles gzip-compressed files transparently.
    Multi-line sequences are concatenated into a single string.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    header = ""
    parts = []
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line
            else:
                parts.append(line)
    return header, "".join(parts)


def write_fasta(path, header, sequence):
    """Write a single-record FASTA file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        if not header.startswith(">"):
            header = f">{header}"
        f.write(f"{header}\n{sequence}\n")


def read_chromosome_lengths(path):
    """Read chromosome lengths from a tab-separated file.

    Expects at least two columns per line: chromosome_name <TAB> length.
    Additional columns (e.g. GenBank accession, RefSeq) are ignored.
    Returns a dict mapping chromosome name to integer length.
    """
    lengths = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    return lengths


def majority_vote(bases, min_freq=1):
    """Return the most frequent valid nucleotide among *bases*, or ``'N'``.

    Only bases in {A, C, G, T} are considered.  Ties are broken in
    alphabetical order (A > C > G > T priority) for reproducibility.
    If no allele reaches *min_freq* occurrences, returns ``'N'``.
    """
    upper = [b.upper() for b in bases]
    counts = {a: upper.count(a) for a in VALID_ALLELES}
    best = max(VALID_ALLELES, key=lambda a: counts[a])
    if counts[best] < min_freq:
        return "N"
    return best


def chrom_id(chrom):
    """Strip a leading ``chr`` prefix, if present.

    Useful for mapping chromosome names between naming conventions,
    e.g. ``chr1`` -> ``1``, ``chrX`` -> ``X``, ``3`` -> ``3``.
    """
    return chrom[3:] if chrom.lower().startswith("chr") else chrom
