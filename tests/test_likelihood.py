"""Tests for ancify.likelihood."""

import numpy as np
import pytest

from ancify.likelihood import (
    BASES,
    JC69,
    K80,
    HKY85,
    GTR,
    build_model,
    call_ancestral_base_likelihood,
    felsenstein_pruning,
)
from ancify.parsimony import parse_newick


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ALL_MODELS = [
    JC69(),
    K80(kappa=2.0),
    HKY85(kappa=2.0, pi=[0.3, 0.2, 0.2, 0.3]),
    GTR(rates=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0], pi=[0.3, 0.2, 0.2, 0.3]),
]

MODEL_IDS = ["JC69", "K80", "HKY85", "GTR"]


# ---------------------------------------------------------------------------
# Rate matrix properties
# ---------------------------------------------------------------------------

class TestRateMatrixProperties:
    """Verify algebraic properties of Q matrices for every model."""

    @pytest.mark.parametrize("model", ALL_MODELS, ids=MODEL_IDS)
    def test_rows_sum_to_zero(self, model):
        Q = model.rate_matrix()
        row_sums = Q.sum(axis=1)
        np.testing.assert_allclose(row_sums, 0.0, atol=1e-12)

    @pytest.mark.parametrize("model", ALL_MODELS, ids=MODEL_IDS)
    def test_off_diagonal_non_negative(self, model):
        Q = model.rate_matrix()
        off_diag = Q.copy()
        np.fill_diagonal(off_diag, 0.0)
        assert np.all(off_diag >= -1e-15)

    @pytest.mark.parametrize("model", ALL_MODELS, ids=MODEL_IDS)
    def test_diagonal_non_positive(self, model):
        Q = model.rate_matrix()
        assert np.all(np.diag(Q) <= 1e-15)

    @pytest.mark.parametrize("model", ALL_MODELS, ids=MODEL_IDS)
    def test_normalized_rate_is_one(self, model):
        Q = model.rate_matrix()
        pi = model.base_frequencies()
        rate = -float(np.dot(pi, np.diag(Q)))
        np.testing.assert_allclose(rate, 1.0, atol=1e-12)

    @pytest.mark.parametrize("model", ALL_MODELS, ids=MODEL_IDS)
    def test_detailed_balance(self, model):
        """π_i Q_ij = π_j Q_ji for all reversible models."""
        Q = model.rate_matrix()
        pi = model.base_frequencies()
        for i in range(4):
            for j in range(4):
                np.testing.assert_allclose(
                    pi[i] * Q[i, j], pi[j] * Q[j, i], atol=1e-12,
                )


# ---------------------------------------------------------------------------
# Transition probability matrix properties
# ---------------------------------------------------------------------------

class TestTransitionProbProperties:

    @pytest.mark.parametrize("model", ALL_MODELS, ids=MODEL_IDS)
    def test_P_zero_is_identity(self, model):
        P = model.transition_probs(0.0)
        np.testing.assert_allclose(P, np.eye(4), atol=1e-12)

    @pytest.mark.parametrize("model", ALL_MODELS, ids=MODEL_IDS)
    def test_rows_sum_to_one(self, model):
        for t in [0.001, 0.01, 0.1, 1.0, 10.0]:
            P = model.transition_probs(t)
            np.testing.assert_allclose(P.sum(axis=1), 1.0, atol=1e-10)

    @pytest.mark.parametrize("model", ALL_MODELS, ids=MODEL_IDS)
    def test_all_entries_non_negative(self, model):
        for t in [0.001, 0.01, 0.1, 1.0]:
            P = model.transition_probs(t)
            assert np.all(P >= -1e-15)

    @pytest.mark.parametrize("model", ALL_MODELS, ids=MODEL_IDS)
    def test_large_t_approaches_equilibrium(self, model):
        P = model.transition_probs(1000.0)
        pi = model.base_frequencies()
        for row in P:
            np.testing.assert_allclose(row, pi, atol=1e-6)

    def test_jc69_closed_form_matches_expm(self):
        jc = JC69()
        for t in [0.01, 0.1, 0.5, 1.0]:
            P_closed = jc.transition_probs(t)
            from scipy.linalg import expm
            P_expm = expm(jc.rate_matrix() * t)
            np.testing.assert_allclose(P_closed, P_expm, atol=1e-10)


# ---------------------------------------------------------------------------
# Model factory
# ---------------------------------------------------------------------------

class TestBuildModel:

    def test_jc69(self):
        m = build_model("JC69")
        assert isinstance(m, JC69)

    def test_k80(self):
        m = build_model("K80", kappa=3.0)
        assert isinstance(m, K80)
        assert m.kappa == 3.0

    def test_hky85(self):
        m = build_model("HKY85", kappa=2.5, base_freqs=[0.3, 0.2, 0.2, 0.3])
        assert isinstance(m, HKY85)

    def test_gtr(self):
        m = build_model("GTR", rates=[1, 2, 1, 1, 2, 1], base_freqs=[0.25] * 4)
        assert isinstance(m, GTR)

    def test_case_insensitive(self):
        assert isinstance(build_model("jc69"), JC69)
        assert isinstance(build_model("hky85"), HKY85)

    def test_unknown_model_raises(self):
        with pytest.raises(ValueError, match="Unknown substitution model"):
            build_model("FAKE_MODEL")


# ---------------------------------------------------------------------------
# Felsenstein pruning
# ---------------------------------------------------------------------------

class TestFelsensteinPruning:

    def test_all_same_base_gives_high_posterior(self):
        tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.15);")
        alleles = {"A": "G", "B": "G", "C": "G"}
        posterior = felsenstein_pruning(tree, alleles, JC69())
        assert posterior["G"] > 0.99

    def test_posteriors_sum_to_one(self):
        tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.15);")
        alleles = {"A": "A", "B": "G", "C": "A"}
        posterior = felsenstein_pruning(tree, alleles, JC69())
        total = sum(posterior.values())
        np.testing.assert_allclose(total, 1.0, atol=1e-10)

    def test_all_missing_gives_uniform(self):
        tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.15);")
        alleles = {"A": "N", "B": "N", "C": "N"}
        posterior = felsenstein_pruning(tree, alleles, JC69())
        for b in BASES:
            np.testing.assert_allclose(posterior[b], 0.25, atol=1e-10)

    def test_single_leaf_data_shifts_posterior(self):
        tree = parse_newick("((A:0.01,B:0.01):0.005,C:0.015);")
        alleles = {"A": "T", "B": "N", "C": "N"}
        posterior = felsenstein_pruning(tree, alleles, JC69())
        assert posterior["T"] > posterior["A"]
        assert posterior["T"] > posterior["C"]
        assert posterior["T"] > posterior["G"]

    def test_distant_outgroup_has_less_influence(self):
        tree = parse_newick("(close:0.01,far:1.0);")
        alleles = {"close": "A", "far": "G"}
        posterior = felsenstein_pruning(tree, alleles, JC69())
        assert posterior["A"] > posterior["G"]

    def test_hky_nonuniform_prior(self):
        tree = parse_newick("(X:0.5,Y:0.5);")
        alleles = {"X": "N", "Y": "N"}
        model = HKY85(kappa=2.0, pi=[0.4, 0.1, 0.1, 0.4])
        posterior = felsenstein_pruning(tree, alleles, model)
        assert posterior["A"] > posterior["C"]
        assert posterior["T"] > posterior["G"]
        np.testing.assert_allclose(posterior["A"], posterior["T"], atol=1e-10)

    def test_works_with_all_models(self):
        tree = parse_newick("((A:0.1,B:0.1):0.05,C:0.15);")
        alleles = {"A": "A", "B": "A", "C": "G"}
        for model in ALL_MODELS:
            posterior = felsenstein_pruning(tree, alleles, model)
            total = sum(posterior.values())
            np.testing.assert_allclose(total, 1.0, atol=1e-10)


# ---------------------------------------------------------------------------
# Confidence encoding
# ---------------------------------------------------------------------------

class TestConfidenceEncoding:
    """Test call_ancestral_base_likelihood confidence-level encoding."""

    TREE = parse_newick(
        "(((bonobo:0.008,chimp:0.008):0.002,gorilla:0.009):0.020,macaque:0.038);"
    )

    def test_all_agree_uppercase(self):
        bases = {"bonobo": "A", "chimp": "A", "gorilla": "A", "macaque": "A"}
        result = call_ancestral_base_likelihood(self.TREE, bases, JC69())
        assert result == "A"
        assert result.isupper()

    def test_all_missing_returns_N(self):
        bases = {"bonobo": "N", "chimp": "N", "gorilla": "N", "macaque": "N"}
        result = call_ancestral_base_likelihood(self.TREE, bases, JC69())
        assert result == "N"

    def test_returns_single_char(self):
        bases = {"bonobo": "G", "chimp": "G", "gorilla": "G", "macaque": "G"}
        result = call_ancestral_base_likelihood(self.TREE, bases, JC69())
        assert len(result) == 1

    def test_result_is_valid_base_or_n(self):
        cases = [
            {"bonobo": "A", "chimp": "A", "gorilla": "A", "macaque": "A"},
            {"bonobo": "N", "chimp": "N", "gorilla": "N", "macaque": "N"},
            {"bonobo": "A", "chimp": "G", "gorilla": "T", "macaque": "C"},
        ]
        for bases in cases:
            result = call_ancestral_base_likelihood(self.TREE, bases, JC69())
            assert result.upper() in ("A", "C", "G", "T", "N")

    def test_low_threshold_gives_lowercase(self):
        bases = {"bonobo": "A", "chimp": "G", "gorilla": "T", "macaque": "C"}
        result = call_ancestral_base_likelihood(
            self.TREE, bases, JC69(), high_threshold=0.99, low_threshold=0.1,
        )
        assert result.islower() or result == "n"

    def test_very_low_threshold_gives_n(self):
        bases = {"bonobo": "A", "chimp": "G", "gorilla": "T", "macaque": "C"}
        result = call_ancestral_base_likelihood(
            self.TREE, bases, JC69(), high_threshold=0.99, low_threshold=0.99,
        )
        assert result == "n"


# ---------------------------------------------------------------------------
# Agreement with parsimony on clear-cut cases
# ---------------------------------------------------------------------------

class TestLikelihoodVsParsimony:
    """On unambiguous cases, likelihood and parsimony should agree on the base."""

    TREE_NWK = "(((bonobo:0.008,chimp:0.008):0.002,gorilla:0.009):0.020,macaque:0.038)"
    TREE = parse_newick(TREE_NWK + ";")

    def test_unanimous_agreement(self):
        from ancify.ancestral import call_ancestral_base_parsimony

        for base in "ACGT":
            alleles = {
                "bonobo": base, "chimp": base,
                "gorilla": base, "macaque": base,
            }
            lik = call_ancestral_base_likelihood(self.TREE, alleles, JC69())
            pars = call_ancestral_base_parsimony(self.TREE, alleles)
            assert lik.upper() == pars.upper() == base

    def test_three_vs_one_agrees_on_majority_base(self):
        from ancify.ancestral import call_ancestral_base_parsimony

        alleles = {"bonobo": "G", "chimp": "G", "gorilla": "A", "macaque": "A"}
        lik = call_ancestral_base_likelihood(self.TREE, alleles, JC69())
        pars = call_ancestral_base_parsimony(self.TREE, alleles)
        assert lik.upper() == pars.upper()
