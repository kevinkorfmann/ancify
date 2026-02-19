"""Tests for vectorized / GPU backend against the scalar reference implementation.

Every test in this module verifies that the vectorized (NumPy / PyTorch)
code paths produce *identical* output to the original per-position
``call_ancestral_base`` and ``majority_vote`` scalar functions.

GPU tests are skipped automatically when no CUDA device is available.
"""

import numpy as np
import pytest

from ancify.ancestral import call_ancestral_base
from ancify.backend import (
    detect_backend,
    encode_sequences,
    get_available_gpus,
    vectorized_ancestral_call,
    vectorized_majority_vote,
)
from ancify.utils import majority_vote

_HAS_GPU = detect_backend() == "gpu"
_DECODE = {0: "A", 1: "C", 2: "G", 3: "T", 4: "N"}
_ENCODE = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}


# ── helpers ──────────────────────────────────────────────────────────────

def _scalar_majority_vote(encoded_matrix, min_freq):
    """Run the scalar majority_vote per-column and return encoded result."""
    n_seqs, L = encoded_matrix.shape
    result = np.empty(L, dtype=np.uint8)
    for col in range(L):
        bases = [_DECODE[encoded_matrix[s, col]] for s in range(n_seqs)]
        result[col] = _ENCODE[majority_vote(bases, min_freq)]
    return result


def _scalar_ancestral_call(inner_seqs, outer_seqs, min_inner, min_outer):
    """Run call_ancestral_base for every position and return the string."""
    length = len(inner_seqs[0])
    out = []
    for i in range(length):
        ib = [s[i].upper() for s in inner_seqs]
        ob = [s[i].upper() for s in outer_seqs]
        out.append(call_ancestral_base(ib, ob, min_inner, min_outer))
    return "".join(out)


def _gpu_device():
    """Return the first CUDA device, or skip the test."""
    import torch
    return torch.device("cuda:0")


# ── encode_sequences ─────────────────────────────────────────────────────

class TestEncodeSequences:

    def test_uppercase(self):
        result = encode_sequences(["ACGT"])
        np.testing.assert_array_equal(result, [[0, 1, 2, 3]])

    def test_lowercase(self):
        result = encode_sequences(["acgt"])
        np.testing.assert_array_equal(result, [[0, 1, 2, 3]])

    def test_mixed_case(self):
        result = encode_sequences(["AcGt"])
        np.testing.assert_array_equal(result, [[0, 1, 2, 3]])

    def test_invalid_chars_become_4(self):
        result = encode_sequences(["ANXM"])
        expected = np.array([[0, 4, 4, 4]], dtype=np.uint8)
        np.testing.assert_array_equal(result, expected)

    def test_multiple_sequences(self):
        result = encode_sequences(["AC", "GT", "NN"])
        expected = np.array([[0, 1], [2, 3], [4, 4]], dtype=np.uint8)
        np.testing.assert_array_equal(result, expected)

    def test_shape(self):
        result = encode_sequences(["ACGT", "TGCA", "AAAA"])
        assert result.shape == (3, 4)
        assert result.dtype == np.uint8


# ── vectorized_majority_vote (CPU) ───────────────────────────────────────

class TestVectorizedMajorityVoteCPU:

    def test_unanimous(self):
        enc = np.array([[0, 1, 2, 3], [0, 1, 2, 3]], dtype=np.uint8)
        result = vectorized_majority_vote(enc, 1)
        np.testing.assert_array_equal(result, [0, 1, 2, 3])

    def test_all_missing(self):
        enc = np.full((3, 5), 4, dtype=np.uint8)
        result = vectorized_majority_vote(enc, 1)
        np.testing.assert_array_equal(result, np.full(5, 4))

    def test_tie_a_beats_c(self):
        enc = np.array([[0], [1]], dtype=np.uint8)
        result = vectorized_majority_vote(enc, 1)
        assert result[0] == 0

    def test_tie_c_beats_g(self):
        enc = np.array([[1], [2]], dtype=np.uint8)
        result = vectorized_majority_vote(enc, 1)
        assert result[0] == 1

    def test_tie_g_beats_t(self):
        enc = np.array([[2], [3]], dtype=np.uint8)
        result = vectorized_majority_vote(enc, 1)
        assert result[0] == 2

    def test_four_way_tie(self):
        enc = np.array([[0], [1], [2], [3]], dtype=np.uint8)
        result = vectorized_majority_vote(enc, 1)
        assert result[0] == 0  # A wins four-way tie

    def test_min_freq_filters(self):
        enc = np.array([[0, 1], [2, 1]], dtype=np.uint8)
        result = vectorized_majority_vote(enc, 2)
        assert result[0] == 4  # A=1,G=1 — neither reaches 2
        assert result[1] == 1  # C=2 — reaches threshold

    def test_min_freq_boundary(self):
        enc = np.array([[0], [0], [1]], dtype=np.uint8)
        assert vectorized_majority_vote(enc, 2)[0] == 0  # A=2 meets min_freq=2
        assert vectorized_majority_vote(enc, 3)[0] == 4  # A=2 fails min_freq=3

    def test_random_matches_scalar(self):
        rng = np.random.RandomState(42)
        enc = rng.randint(0, 5, size=(5, 2000)).astype(np.uint8)
        for min_freq in [1, 2, 3]:
            vec = vectorized_majority_vote(enc, min_freq)
            scalar = _scalar_majority_vote(enc, min_freq)
            np.testing.assert_array_equal(
                vec, scalar, err_msg=f"min_freq={min_freq}"
            )

    def test_random_large_matches_scalar(self):
        rng = np.random.RandomState(7)
        enc = rng.randint(0, 5, size=(10, 5000)).astype(np.uint8)
        for min_freq in [1, 4, 8]:
            vec = vectorized_majority_vote(enc, min_freq)
            scalar = _scalar_majority_vote(enc, min_freq)
            np.testing.assert_array_equal(
                vec, scalar, err_msg=f"min_freq={min_freq}"
            )


# ── vectorized_ancestral_call (CPU) ──────────────────────────────────────

class TestVectorizedAncestralCallCPU:

    def test_all_agree_uppercase(self):
        inner = ["ACGT", "ACGT", "ACGT"]
        outer = ["ACGT"]
        assert vectorized_ancestral_call(inner, outer) == "ACGT"

    def test_inner_only_lowercase(self):
        inner = ["ACGT", "ACGT"]
        outer = ["NNNN"]
        assert vectorized_ancestral_call(inner, outer) == "acgt"

    def test_outer_only_lowercase(self):
        inner = ["NNNN"]
        outer = ["ACGT"]
        assert vectorized_ancestral_call(inner, outer) == "acgt"

    def test_both_missing_gives_N(self):
        inner = ["NNNN"]
        outer = ["NNNN"]
        assert vectorized_ancestral_call(inner, outer) == "NNNN"

    def test_disagree_gives_n(self):
        inner = ["AAAA"]
        outer = ["CCCC"]
        assert vectorized_ancestral_call(inner, outer) == "nnnn"

    def test_mixed_confidence(self):
        inner = ["ACNN"]
        outer = ["ANNG"]
        result = vectorized_ancestral_call(inner, outer)
        assert result[0] == "A"  # agree -> uppercase
        assert result[1] == "c"  # inner only -> lowercase
        assert result[2] == "N"  # both missing
        assert result[3] == "g"  # outer only -> lowercase

    def test_inner_disagree_with_outer(self):
        inner = ["ACGT"]
        outer = ["TGCA"]
        result = vectorized_ancestral_call(inner, outer)
        assert result == "nnnn"

    def test_min_inner_freq_filters(self):
        inner = ["A", "G", "C"]
        outer = ["A"]
        result = vectorized_ancestral_call(
            inner, outer, min_inner_freq=2, min_outer_freq=1
        )
        assert result == "a"

    def test_min_inner_freq_passes(self):
        inner = ["A", "A", "G"]
        outer = ["A"]
        result = vectorized_ancestral_call(
            inner, outer, min_inner_freq=2, min_outer_freq=1
        )
        assert result == "A"

    def test_min_outer_freq_filters(self):
        inner = ["A", "A"]
        outer = ["A"]
        result = vectorized_ancestral_call(
            inner, outer, min_inner_freq=1, min_outer_freq=2
        )
        assert result == "a"

    def test_random_matches_scalar(self):
        rng = np.random.RandomState(99)
        bases = list("ACGTN")
        length = 500
        for n_inner, n_outer in [(3, 1), (3, 2), (1, 1), (5, 3)]:
            inner = [
                "".join(rng.choice(bases, size=length)) for _ in range(n_inner)
            ]
            outer = [
                "".join(rng.choice(bases, size=length)) for _ in range(n_outer)
            ]
            for min_i in [1, 2]:
                for min_o in [1, 2]:
                    vec = vectorized_ancestral_call(
                        inner, outer, min_i, min_o
                    )
                    scalar = _scalar_ancestral_call(
                        inner, outer, min_i, min_o
                    )
                    assert vec == scalar, (
                        f"Mismatch: n_inner={n_inner} n_outer={n_outer} "
                        f"min_i={min_i} min_o={min_o}"
                    )

    def test_all_bases_high_confidence(self):
        for base in "ACGT":
            inner = [base * 10] * 3
            outer = [base * 10]
            result = vectorized_ancestral_call(inner, outer)
            assert result == base * 10

    def test_all_bases_low_confidence_inner_only(self):
        for base in "ACGT":
            inner = [base * 5] * 2
            outer = ["N" * 5]
            result = vectorized_ancestral_call(inner, outer)
            assert result == base.lower() * 5

    def test_all_bases_low_confidence_outer_only(self):
        for base in "ACGT":
            inner = ["N" * 5]
            outer = [base * 5]
            result = vectorized_ancestral_call(inner, outer)
            assert result == base.lower() * 5


# ── GPU tests ────────────────────────────────────────────────────────────

@pytest.mark.skipif(not _HAS_GPU, reason="No CUDA GPU available")
class TestVectorizedMajorityVoteGPU:

    def test_unanimous(self):
        device = _gpu_device()
        enc = np.array([[0, 1, 2, 3], [0, 1, 2, 3]], dtype=np.uint8)
        result = vectorized_majority_vote(enc, 1, device=device)
        np.testing.assert_array_equal(result, [0, 1, 2, 3])

    def test_tie_breaking(self):
        device = _gpu_device()
        enc = np.array([[0], [1], [2], [3]], dtype=np.uint8)
        result = vectorized_majority_vote(enc, 1, device=device)
        assert result[0] == 0  # A wins four-way tie

    def test_matches_cpu(self):
        device = _gpu_device()
        rng = np.random.RandomState(123)
        enc = rng.randint(0, 5, size=(6, 3000)).astype(np.uint8)
        for min_freq in [1, 2, 4]:
            cpu = vectorized_majority_vote(enc, min_freq, device=None)
            gpu = vectorized_majority_vote(enc, min_freq, device=device)
            np.testing.assert_array_equal(
                cpu, gpu, err_msg=f"min_freq={min_freq}"
            )

    def test_matches_scalar(self):
        device = _gpu_device()
        rng = np.random.RandomState(200)
        enc = rng.randint(0, 5, size=(5, 2000)).astype(np.uint8)
        for min_freq in [1, 3]:
            gpu = vectorized_majority_vote(enc, min_freq, device=device)
            scalar = _scalar_majority_vote(enc, min_freq)
            np.testing.assert_array_equal(
                gpu, scalar, err_msg=f"min_freq={min_freq}"
            )


@pytest.mark.skipif(not _HAS_GPU, reason="No CUDA GPU available")
class TestVectorizedAncestralCallGPU:

    def test_all_agree_uppercase(self):
        device = _gpu_device()
        inner = ["ACGT", "ACGT"]
        outer = ["ACGT"]
        result = vectorized_ancestral_call(inner, outer, device=device)
        assert result == "ACGT"

    def test_mixed_confidence(self):
        device = _gpu_device()
        inner = ["ACNN"]
        outer = ["ANNG"]
        result = vectorized_ancestral_call(inner, outer, device=device)
        assert result[0] == "A"
        assert result[1] == "c"
        assert result[2] == "N"
        assert result[3] == "g"

    def test_matches_cpu(self):
        device = _gpu_device()
        rng = np.random.RandomState(77)
        bases = list("ACGTN")
        length = 2000
        inner = ["".join(rng.choice(bases, size=length)) for _ in range(3)]
        outer = ["".join(rng.choice(bases, size=length)) for _ in range(2)]
        cpu = vectorized_ancestral_call(inner, outer, device=None)
        gpu = vectorized_ancestral_call(inner, outer, device=device)
        assert cpu == gpu

    def test_matches_scalar(self):
        device = _gpu_device()
        rng = np.random.RandomState(55)
        bases = list("ACGTN")
        length = 1000
        inner = ["".join(rng.choice(bases, size=length)) for _ in range(4)]
        outer = ["".join(rng.choice(bases, size=length)) for _ in range(1)]
        gpu = vectorized_ancestral_call(inner, outer, device=device)
        scalar = _scalar_ancestral_call(inner, outer, 1, 1)
        assert gpu == scalar

    def test_matches_scalar_with_thresholds(self):
        device = _gpu_device()
        rng = np.random.RandomState(33)
        bases = list("ACGTN")
        length = 800
        inner = ["".join(rng.choice(bases, size=length)) for _ in range(5)]
        outer = ["".join(rng.choice(bases, size=length)) for _ in range(3)]
        for min_i in [1, 2, 3]:
            for min_o in [1, 2]:
                gpu = vectorized_ancestral_call(
                    inner, outer, min_i, min_o, device=device
                )
                scalar = _scalar_ancestral_call(inner, outer, min_i, min_o)
                assert gpu == scalar, f"min_i={min_i}, min_o={min_o}"


# ── Backend detection (smoke tests) ─────────────────────────────────────

class TestBackendDetection:

    def test_detect_returns_str(self):
        result = detect_backend()
        assert result in ("cpu", "gpu")

    def test_get_available_gpus_returns_list(self):
        result = get_available_gpus()
        assert isinstance(result, list)
        for item in result:
            assert isinstance(item, int)
