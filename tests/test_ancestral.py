"""Tests for ancify.ancestral."""

import pytest

from ancify.ancestral import call_ancestral_base, call_ancestral_base_parsimony
from ancify.parsimony import parse_newick


class TestCallAncestralBase:
    """Test the core ancestral calling function with various scenarios."""

    # ── High confidence (inner == outer) ──────────────────────────────────

    def test_all_agree(self):
        assert call_ancestral_base(["A", "A", "A"], ["A"]) == "A"

    def test_inner_majority_agrees_with_outer(self):
        assert call_ancestral_base(["A", "A", "G"], ["A"]) == "A"

    def test_all_species_same_base(self):
        for base in "ACGT":
            assert call_ancestral_base([base, base], [base]) == base

    def test_high_confidence_is_uppercase(self):
        result = call_ancestral_base(["T", "T", "T"], ["T"])
        assert result == "T"
        assert result.isupper()

    # ── Low confidence (one tier missing) ─────────────────────────────────

    def test_inner_present_outer_missing(self):
        result = call_ancestral_base(["A", "A", "A"], ["N"])
        assert result == "a"
        assert result.islower()

    def test_outer_present_inner_missing(self):
        result = call_ancestral_base(["N", "N", "N"], ["G"])
        assert result == "g"
        assert result.islower()

    def test_inner_present_outer_all_gaps(self):
        result = call_ancestral_base(["C", "C"], ["-", "-"])
        assert result == "c"

    def test_outer_present_inner_all_gaps(self):
        result = call_ancestral_base(["-", "-", "-"], ["T"])
        assert result == "t"

    # ── Missing (both tiers lack data) ────────────────────────────────────

    def test_both_missing(self):
        assert call_ancestral_base(["N", "N"], ["N"]) == "N"

    def test_both_gaps(self):
        assert call_ancestral_base(["-"], ["-"]) == "N"

    def test_all_invalid(self):
        assert call_ancestral_base(["N", "-", "N"], ["N", "-"]) == "N"

    # ── Disagreement (inner != outer, both non-N) ─────────────────────────

    def test_disagree(self):
        result = call_ancestral_base(["A", "A", "A"], ["T"])
        assert result == "n"

    def test_disagree_different_bases(self):
        assert call_ancestral_base(["C", "C"], ["G"]) == "n"

    # ── min_freq thresholds ───────────────────────────────────────────────

    def test_min_inner_freq_too_high(self):
        result = call_ancestral_base(
            ["A", "G", "C"], ["A"], min_inner_freq=2
        )
        assert result == "a"  # inner fails threshold -> low confidence from outer

    def test_min_inner_freq_met(self):
        result = call_ancestral_base(
            ["A", "A", "G"], ["A"], min_inner_freq=2
        )
        assert result == "A"

    def test_min_outer_freq_too_high(self):
        result = call_ancestral_base(
            ["A", "A"], ["A"], min_outer_freq=2
        )
        assert result == "a"  # outer fails threshold -> low confidence from inner

    def test_min_outer_freq_met(self):
        result = call_ancestral_base(
            ["A", "A"], ["A", "A"], min_outer_freq=2
        )
        assert result == "A"

    # ── Edge cases ────────────────────────────────────────────────────────

    def test_single_inner_single_outer(self):
        assert call_ancestral_base(["G"], ["G"]) == "G"

    def test_many_inner_species(self):
        inner = ["A"] * 10 + ["G"] * 3
        assert call_ancestral_base(inner, ["A"]) == "A"

    def test_many_outer_species(self):
        outer = ["T"] * 5
        assert call_ancestral_base(["T", "T"], outer) == "T"

    def test_inner_tie_with_outer_matching_one(self):
        # Inner tie between A and C (1 each), majority_vote picks A (alphabetical)
        result = call_ancestral_base(["A", "C"], ["A"])
        assert result == "A"

    def test_inner_tie_with_outer_matching_other(self):
        result = call_ancestral_base(["A", "C"], ["C"])
        # Inner majority is A (tie-break), outer is C -> disagreement
        assert result == "n"

    def test_lowercase_input_handled(self):
        result = call_ancestral_base(["a", "a", "a"], ["a"])
        assert result == "A"


class TestCallAncestralBaseParsimony:
    """Test the parsimony-based ancestral calling function."""

    TREE = parse_newick("(((bonobo,chimp),gorilla),macaque);")

    def test_all_agree_uppercase(self):
        bases = {"bonobo": "A", "chimp": "A", "gorilla": "A", "macaque": "A"}
        result = call_ancestral_base_parsimony(self.TREE, bases)
        assert result == "A"
        assert result.isupper()

    def test_ambiguous_lowercase(self):
        bases = {"bonobo": "A", "chimp": "A", "gorilla": "A", "macaque": "G"}
        result = call_ancestral_base_parsimony(self.TREE, bases)
        assert result.islower()
        assert result.upper() in "ACGT"

    def test_all_missing_returns_N(self):
        bases = {"bonobo": "N", "chimp": "N", "gorilla": "N", "macaque": "N"}
        result = call_ancestral_base_parsimony(self.TREE, bases)
        assert result == "N"

    def test_some_missing_still_resolves(self):
        bases = {"bonobo": "T", "chimp": "T", "gorilla": "N", "macaque": "T"}
        result = call_ancestral_base_parsimony(self.TREE, bases)
        assert result == "T"

    def test_simple_two_taxa(self):
        tree = parse_newick("(X,Y);")
        bases = {"X": "C", "Y": "C"}
        result = call_ancestral_base_parsimony(tree, bases)
        assert result == "C"

    def test_two_taxa_disagree_lowercase(self):
        tree = parse_newick("(X,Y);")
        bases = {"X": "A", "Y": "T"}
        result = call_ancestral_base_parsimony(tree, bases)
        assert result.islower()
        assert result.upper() in "ACGT"

    def test_returns_single_char(self):
        bases = {"bonobo": "G", "chimp": "G", "gorilla": "G", "macaque": "G"}
        result = call_ancestral_base_parsimony(self.TREE, bases)
        assert len(result) == 1

    def test_result_is_valid_or_N(self):
        cases = [
            {"bonobo": "A", "chimp": "A", "gorilla": "A", "macaque": "A"},
            {"bonobo": "N", "chimp": "N", "gorilla": "N", "macaque": "N"},
            {"bonobo": "A", "chimp": "G", "gorilla": "T", "macaque": "C"},
        ]
        for bases in cases:
            result = call_ancestral_base_parsimony(self.TREE, bases)
            assert result.upper() in ("A", "C", "G", "T", "N")

    def test_missing_species_treated_as_N(self):
        bases = {"bonobo": "A", "chimp": "A"}
        result = call_ancestral_base_parsimony(self.TREE, bases)
        assert result.upper() in "ACGT"

    def test_lowercase_input_handled(self):
        bases = {"bonobo": "c", "chimp": "c", "gorilla": "c", "macaque": "c"}
        result = call_ancestral_base_parsimony(self.TREE, bases)
        assert result == "C"


class TestCallAncestralBaseReturnTypes:
    """Verify return type properties."""

    def test_always_returns_single_char(self):
        cases = [
            (["A", "A"], ["A"]),
            (["N", "N"], ["N"]),
            (["A", "A"], ["G"]),
            (["N"], ["T"]),
            (["C"], ["N"]),
        ]
        for inner, outer in cases:
            result = call_ancestral_base(inner, outer)
            assert len(result) == 1

    def test_result_is_valid_or_n(self):
        cases = [
            (["A", "A"], ["A"]),
            (["N", "N"], ["N"]),
            (["A", "A"], ["G"]),
            (["N"], ["T"]),
        ]
        for inner, outer in cases:
            result = call_ancestral_base(inner, outer)
            assert result.upper() in ("A", "C", "G", "T", "N")
