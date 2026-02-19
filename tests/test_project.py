"""Tests for ancify.project."""

import textwrap

import pytest

from ancify.project import project_alignment


class TestProjectAlignment:
    def test_basic_projection(self, minimal_axt):
        """Bases from the AXT outgroup line land at the correct focal positions."""
        seq = project_alignment(str(minimal_axt), chrom_length=30, target_chrom="chr1")
        assert len(seq) == 30
        # Positions 0-9 (0-based) from alignment block 0 (chr1 1-10, 1-based)
        assert seq[0] == "A"
        assert seq[1] == "C"
        assert seq[8] == "A"  # outgroup had 'A' where focal had 'N'
        assert seq[9] == "A"  # outgroup had 'A' where focal had 'N'
        # Positions 14-18 (0-based) from alignment block 1 (chr1 15-19)
        assert seq[14] == "G"
        assert seq[15] == "G"

    def test_unaligned_positions_are_N(self, minimal_axt):
        """Positions not covered by any alignment block should be N."""
        seq = project_alignment(str(minimal_axt), chrom_length=30, target_chrom="chr1")
        assert seq[10] == "N"
        assert seq[11] == "N"
        assert seq[25] == "N"
        assert seq[29] == "N"

    def test_wrong_chromosome_all_N(self, minimal_axt):
        """Querying a chromosome not present in the AXT gives all N."""
        seq = project_alignment(str(minimal_axt), chrom_length=20, target_chrom="chr99")
        assert seq == "N" * 20

    def test_gaps_in_focal_are_skipped(self, axt_with_gaps):
        """Gaps in the focal sequence (insertions in outgroup) don't advance position."""
        seq = project_alignment(str(axt_with_gaps), chrom_length=20, target_chrom="chr1")
        assert len(seq) == 20
        # Focal: AC-GTA-CGTN  -> non-gap positions: A C G T A C G T N  (9 focal positions)
        # Query: ACAGTANCGTN
        # Position 0: A->A, 1: C->C, 2(gap skip), 2: G->G, 3: T->T, 4: A->A
        # 5(gap skip), 5: C->C, 6: G->G, 7: T->T, 8: N->N
        assert seq[0] == "A"
        assert seq[1] == "C"
        assert seq[2] == "G"  # focal 'G' aligned to outgroup 'G' (after gap skip)
        assert seq[3] == "T"
        assert seq[4] == "A"  # focal 'A' aligned to outgroup 'N' (wait, let me recheck)

    def test_gzipped_axt(self, tmp_path):
        """Gzipped AXT files are read correctly."""
        import gzip
        content = textwrap.dedent("""\
            0 chr1 1 5 chrQ 1 5 + 100
            ACGTA
            TTTTT

        """)
        p = tmp_path / "test.axt.gz"
        with gzip.open(p, "wt") as f:
            f.write(content)
        seq = project_alignment(str(p), chrom_length=10, target_chrom="chr1")
        assert seq[0] == "T"
        assert seq[4] == "T"
        assert seq[5] == "N"

    def test_empty_axt(self, tmp_path):
        """An empty AXT file produces all N."""
        p = tmp_path / "empty.axt"
        p.write_text("")
        seq = project_alignment(str(p), chrom_length=10, target_chrom="chr1")
        assert seq == "N" * 10

    def test_preserves_outgroup_case(self, tmp_path):
        """Outgroup base case (upper/lower) is preserved in projection."""
        content = textwrap.dedent("""\
            0 chr1 1 4 chrQ 1 4 + 100
            ACGT
            AcGt

        """)
        p = tmp_path / "case.axt"
        p.write_text(content)
        seq = project_alignment(str(p), chrom_length=10, target_chrom="chr1")
        assert seq[0] == "A"
        assert seq[1] == "c"
        assert seq[2] == "G"
        assert seq[3] == "t"

    def test_multiple_blocks_same_chrom(self, tmp_path):
        """Multiple alignment blocks for the same chromosome are merged."""
        content = textwrap.dedent("""\
            0 chr1 1 3 chrQ 1 3 + 100
            AAA
            TTT

            1 chr1 8 10 chrQ 50 52 + 100
            GGG
            CCC

        """)
        p = tmp_path / "multi.axt"
        p.write_text(content)
        seq = project_alignment(str(p), chrom_length=15, target_chrom="chr1")
        assert seq[0:3] == "TTT"
        assert seq[3:7] == "NNNN"
        assert seq[7:10] == "CCC"
        assert seq[10:] == "NNNNN"

    def test_comment_lines_ignored(self, tmp_path):
        """Lines starting with # are skipped."""
        content = textwrap.dedent("""\
            # This is a comment
            0 chr1 1 3 chrQ 1 3 + 100
            ACG
            TTT

        """)
        p = tmp_path / "comments.axt"
        p.write_text(content)
        seq = project_alignment(str(p), chrom_length=5, target_chrom="chr1")
        assert seq[0:3] == "TTT"
