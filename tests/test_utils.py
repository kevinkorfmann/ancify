"""Tests for ancify.utils."""

import gzip
from pathlib import Path

import pytest

from ancify.utils import (
    VALID_ALLELES,
    chrom_id,
    majority_vote,
    read_chromosome_lengths,
    read_fasta,
    write_fasta,
)


# ── read_fasta ────────────────────────────────────────────────────────────────

class TestReadFasta:
    def test_simple(self, simple_fasta):
        header, seq = read_fasta(simple_fasta)
        assert header == ">chr1"
        assert seq == "ACGTACGTNN"

    def test_multiline(self, multiline_fasta):
        header, seq = read_fasta(multiline_fasta)
        assert header == ">chrX"
        assert seq == "ACGTNNNNTTTT"

    def test_gzipped(self, tmp_path):
        p = tmp_path / "test.fa.gz"
        with gzip.open(p, "wt") as f:
            f.write(">chrZ\nACGT\nTTTT\n")
        header, seq = read_fasta(p)
        assert header == ">chrZ"
        assert seq == "ACGTTTTT"

    def test_empty_sequence(self, tmp_path):
        p = tmp_path / "empty.fa"
        p.write_text(">chr1\n")
        header, seq = read_fasta(p)
        assert header == ">chr1"
        assert seq == ""

    def test_no_header(self, tmp_path):
        p = tmp_path / "noheader.fa"
        p.write_text("ACGTACGT\n")
        header, seq = read_fasta(p)
        assert header == ""
        assert seq == "ACGTACGT"


# ── write_fasta ───────────────────────────────────────────────────────────────

class TestWriteFasta:
    def test_roundtrip(self, tmp_path):
        p = tmp_path / "out.fa"
        write_fasta(p, ">myseq", "ACGTNNNN")
        header, seq = read_fasta(p)
        assert header == ">myseq"
        assert seq == "ACGTNNNN"

    def test_adds_chevron(self, tmp_path):
        p = tmp_path / "out.fa"
        write_fasta(p, "myseq", "ACGT")
        header, seq = read_fasta(p)
        assert header == ">myseq"

    def test_creates_parent_dirs(self, tmp_path):
        p = tmp_path / "a" / "b" / "c" / "out.fa"
        write_fasta(p, ">deep", "TTTT")
        assert p.exists()
        _, seq = read_fasta(p)
        assert seq == "TTTT"


# ── read_chromosome_lengths ──────────────────────────────────────────────────

class TestReadChromosomeLengths:
    def test_basic(self, chrom_lengths_file):
        lengths = read_chromosome_lengths(chrom_lengths_file)
        assert lengths == {"chr1": 100, "chr2": 80, "chrX": 60}

    def test_two_column(self, tmp_path):
        p = tmp_path / "lens.txt"
        p.write_text("2L\t23513712\n2R\t25286936\n")
        lengths = read_chromosome_lengths(p)
        assert lengths == {"2L": 23513712, "2R": 25286936}

    def test_empty_file(self, tmp_path):
        p = tmp_path / "empty.txt"
        p.write_text("")
        lengths = read_chromosome_lengths(p)
        assert lengths == {}

    def test_blank_lines_skipped(self, tmp_path):
        p = tmp_path / "blanks.txt"
        p.write_text("chr1\t100\n\nchr2\t200\n\n")
        lengths = read_chromosome_lengths(p)
        assert lengths == {"chr1": 100, "chr2": 200}


# ── majority_vote ────────────────────────────────────────────────────────────

class TestMajorityVote:
    def test_unanimous(self):
        assert majority_vote(["A", "A", "A"]) == "A"

    def test_majority(self):
        assert majority_vote(["A", "A", "G"]) == "A"

    def test_all_different(self):
        result = majority_vote(["A", "C", "G"])
        assert result in VALID_ALLELES

    def test_tie_breaks_alphabetically(self):
        assert majority_vote(["A", "C"]) == "A"
        assert majority_vote(["C", "G"]) == "C"
        assert majority_vote(["G", "T"]) == "G"

    def test_all_missing(self):
        assert majority_vote(["N", "N", "N"]) == "N"

    def test_single_species(self):
        assert majority_vote(["T"]) == "T"

    def test_min_freq_threshold(self):
        assert majority_vote(["A", "G", "C"], min_freq=2) == "N"
        assert majority_vote(["A", "A", "G"], min_freq=2) == "A"

    def test_min_freq_equals_count(self):
        assert majority_vote(["A", "A", "A"], min_freq=3) == "A"
        assert majority_vote(["A", "A", "G"], min_freq=3) == "N"

    def test_ignores_N_and_gaps(self):
        assert majority_vote(["A", "N", "-"]) == "A"

    def test_lowercase_input(self):
        assert majority_vote(["a", "a", "g"]) == "A"

    def test_empty_list(self):
        assert majority_vote([]) == "N"

    def test_only_gaps(self):
        assert majority_vote(["-", "-", "-"]) == "N"


# ── chrom_id ─────────────────────────────────────────────────────────────────

class TestChromId:
    def test_strip_chr(self):
        assert chrom_id("chr1") == "1"
        assert chrom_id("chrX") == "X"
        assert chrom_id("chr22") == "22"

    def test_no_prefix(self):
        assert chrom_id("2L") == "2L"
        assert chrom_id("X") == "X"
        assert chrom_id("MT") == "MT"

    def test_case_insensitive_prefix(self):
        assert chrom_id("Chr1") == "1"
        assert chrom_id("CHR1") == "1"

    def test_empty(self):
        assert chrom_id("chr") == ""
