"""Shared fixtures for ancify tests."""

import textwrap
from pathlib import Path

import pytest


@pytest.fixture
def tmp_dir(tmp_path):
    """Return a temporary directory as a Path."""
    return tmp_path


@pytest.fixture
def simple_fasta(tmp_path):
    """Create a simple single-record FASTA file and return its path."""
    p = tmp_path / "test.fa"
    p.write_text(">chr1\nACGTACGTNN\n")
    return p


@pytest.fixture
def multiline_fasta(tmp_path):
    """Create a multi-line FASTA file and return its path."""
    p = tmp_path / "multi.fa"
    p.write_text(">chrX\nACGT\nNNNN\nTTTT\n")
    return p


@pytest.fixture
def chrom_lengths_file(tmp_path):
    """Create a chromosome-lengths file and return its path."""
    p = tmp_path / "chromoLens.txt"
    p.write_text("chr1\t100\tCM000001\tNC_000001\nchr2\t80\tCM000002\tNC_000002\nchrX\t60\n")
    return p


@pytest.fixture
def minimal_axt(tmp_path):
    """Create a minimal AXT alignment file.

    This represents a tiny alignment: focal chr1 positions 1-10 (1-based)
    aligned to an outgroup with known bases.
    """
    content = textwrap.dedent("""\
        0 chr1 1 10 chrUn 500 509 + 1000
        ACGTACGTNN
        ACGTACGTAA

        1 chr1 15 19 chrUn 600 604 + 500
        TTTTT
        GGGGG

    """)
    p = tmp_path / "test.axt"
    p.write_text(content)
    return p


@pytest.fixture
def axt_with_gaps(tmp_path):
    """AXT alignment with gaps in the focal sequence (outgroup insertions)."""
    content = textwrap.dedent("""\
        0 chr1 1 7 chrUn 100 109 + 800
        AC-GTA-CGTN
        ACAGTANCGTN

    """)
    p = tmp_path / "test_gaps.axt"
    p.write_text(content)
    return p


@pytest.fixture
def projected_fastas(tmp_path):
    """Create a set of projected FASTA files for 3 inner + 1 outer species.

    Sequence length = 20.  Inner species: sp1, sp2, sp3.  Outer: sp4.
    """
    seqs = {
        "sp1": "ACGTACGTNNNNNNNNNNN",
        "sp2": "ACGTACGTNNNNNNNNNNN",
        "sp3": "ACGTACGTNNNNNNNNNNN",
        "sp4": "ACGTNNNNNNNNNNNNNN",
    }
    proj = tmp_path / "projected"
    for name, seq in seqs.items():
        d = proj / name
        d.mkdir(parents=True)
        (d / "chr1.fa").write_text(f">chr1\n{seq}\n")
    return tmp_path, seqs


@pytest.fixture
def full_pipeline_dir(tmp_path):
    """Set up a complete mini-pipeline directory for integration testing.

    Creates AXT files, chromosome lengths, and config for a tiny genome.
    """
    chrom_len = 30

    lengths = tmp_path / "chromoLens.txt"
    lengths.write_text(f"chr1\t{chrom_len}\n")

    axt_template = textwrap.dedent("""\
        0 chr1 1 10 chrQ 1 10 + 100
        ACGTACGTAC
        {seq}

    """)

    species = {
        "inner1": {"seq": "ACGTACGTAC", "file": "focal.inner1.axt"},
        "inner2": {"seq": "ACGTACGTAC", "file": "focal.inner2.axt"},
        "outer1": {"seq": "ACGTACGTAC", "file": "focal.outer1.axt"},
    }
    for name, info in species.items():
        p = tmp_path / info["file"]
        p.write_text(axt_template.format(seq=info["seq"]))

    config_text = textwrap.dedent(f"""\
        focal_species: test_species
        chromosome_lengths: {lengths}
        chromosomes:
          - chr1
        outgroups:
          inner:
            - name: inner1
              alignment: {tmp_path / 'focal.inner1.axt'}
            - name: inner2
              alignment: {tmp_path / 'focal.inner2.axt'}
          outer:
            - name: outer1
              alignment: {tmp_path / 'focal.outer1.axt'}
        work_dir: {tmp_path}
        output_dir: {tmp_path / 'ancestral_calls'}
        min_inner_freq: 1
        min_outer_freq: 1
        num_cpus: 1
    """)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(config_text)

    return tmp_path, config_path
