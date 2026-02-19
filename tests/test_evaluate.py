"""Tests for ancify.evaluate."""

import pytest
import numpy as np

from ancify.evaluate import (
    compute_coverage_stats,
    compare_to_reference,
    _format_pattern,
)


class TestComputeCoverageStats:
    def test_all_high_confidence(self):
        stats = compute_coverage_stats("ACGTACGT")
        assert stats["total_positions"] == 8
        assert stats["high_confidence"] == 8
        assert stats["low_confidence"] == 0
        assert stats["missing"] == 0
        assert stats["prop_nonmissing"] == pytest.approx(1.0)
        assert stats["prop_high_confidence"] == pytest.approx(1.0)

    def test_all_missing(self):
        stats = compute_coverage_stats("NNNNNNNN")
        assert stats["high_confidence"] == 0
        assert stats["low_confidence"] == 0
        assert stats["missing"] == 8
        assert stats["prop_nonmissing"] == pytest.approx(0.0)

    def test_mixed_confidence(self):
        # 4 uppercase (high), 2 lowercase (low), 2 N (missing)
        stats = compute_coverage_stats("ACGTacNN")
        assert stats["total_positions"] == 8
        assert stats["high_confidence"] == 4
        assert stats["low_confidence"] == 2
        assert stats["missing"] == 2
        assert stats["prop_nonmissing"] == pytest.approx(6 / 8)
        assert stats["prop_high_confidence"] == pytest.approx(4 / 8)

    def test_lowercase_n_is_missing(self):
        stats = compute_coverage_stats("ACnn")
        assert stats["high_confidence"] == 2
        assert stats["low_confidence"] == 0
        assert stats["missing"] == 2

    def test_single_position(self):
        stats = compute_coverage_stats("A")
        assert stats["total_positions"] == 1
        assert stats["high_confidence"] == 1

    def test_all_low_confidence(self):
        stats = compute_coverage_stats("acgt")
        assert stats["high_confidence"] == 0
        assert stats["low_confidence"] == 4
        assert stats["missing"] == 0


class TestCompareToReference:
    def test_identical_sequences(self):
        stats = compare_to_reference("ACGTACGT", "ACGTACGT")
        assert stats["agreement_rate"] == pytest.approx(1.0)
        assert stats["disagreement_rate"] == pytest.approx(0.0)

    def test_complete_disagreement(self):
        stats = compare_to_reference("AAAA", "TTTT")
        assert stats["agreement_rate"] == pytest.approx(0.0)
        assert stats["disagreement_rate"] == pytest.approx(1.0)

    def test_partial_agreement(self):
        stats = compare_to_reference("AAGT", "AACT")
        # 3 of 4 agree
        assert stats["agreement_rate"] == pytest.approx(0.75)

    def test_missing_in_test(self):
        stats = compare_to_reference("AANN", "AAGT")
        assert stats["test_nonmissing"] == pytest.approx(0.5)
        assert stats["ref_nonmissing"] == pytest.approx(1.0)
        assert stats["both_nonmissing"] == 2
        assert stats["agreement_rate"] == pytest.approx(1.0)

    def test_missing_in_reference(self):
        stats = compare_to_reference("ACGT", "ACNN")
        assert stats["both_nonmissing"] == 2
        assert stats["agreement_rate"] == pytest.approx(1.0)

    def test_case_insensitive_comparison(self):
        stats = compare_to_reference("acgt", "ACGT")
        assert stats["agreement_rate"] == pytest.approx(1.0)

    def test_all_missing(self):
        stats = compare_to_reference("NNNN", "NNNN")
        assert stats["both_nonmissing"] == 0
        assert np.isnan(stats["agreement_rate"])

    def test_with_positions_subset(self):
        test = "ACGTACGT"
        ref = "ACGTTTTT"
        positions = np.array([0, 1, 2, 3])
        stats = compare_to_reference(test, ref, positions=positions)
        assert stats["agreement_rate"] == pytest.approx(1.0)


class TestFormatPattern:
    def test_chrom_placeholder(self):
        assert _format_pattern("file_{chrom}.fa", "chr1") == "file_chr1.fa"

    def test_chrom_id_placeholder(self):
        assert _format_pattern("file_{chrom_id}.fa", "chr1") == "file_1.fa"

    def test_both_placeholders(self):
        result = _format_pattern("{chrom}_{chrom_id}.fa", "chrX")
        assert result == "chrX_X.fa"

    def test_no_chr_prefix(self):
        assert _format_pattern("file_{chrom_id}.fa", "2L") == "file_2L.fa"
        assert _format_pattern("file_{chrom}.fa", "2L") == "file_2L.fa"

    def test_no_placeholder(self):
        assert _format_pattern("static.fa", "chr1") == "static.fa"
