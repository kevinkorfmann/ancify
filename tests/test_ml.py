"""Tests for ancify.ml -- ML-based ancestral allele calling."""

import numpy as np
import pytest

from ancify.ml import (
    _cpg_flag,
    _gc_content,
    _labels_from_reference,
    _labels_from_voting,
    _stratified_subsample,
    encode_confidence,
    extract_features,
    generate_training_data,
    predict,
)

try:
    import lightgbm  # noqa: F401
    HAS_LIGHTGBM = True
except ImportError:
    HAS_LIGHTGBM = False

requires_lightgbm = pytest.mark.skipif(
    not HAS_LIGHTGBM, reason="lightgbm not installed"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_seqs(inner, outer):
    """Pad short sequences to equal length if needed."""
    return inner, outer


# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------

class TestExtractFeatures:
    """Verify shape, dtype, and semantics of the feature matrix."""

    def test_shape_basic(self):
        inner = ["ACGT", "ACGT"]
        outer = ["ACGT"]
        feat = extract_features(inner, outer)
        K = len(inner) + len(outer)
        expected_cols = K + 9  # K outgroup bases + 9 derived features
        assert feat.shape == (4, expected_cols)
        assert feat.dtype == np.float32

    def test_shape_longer_sequence(self):
        L = 200
        inner = ["A" * L, "C" * L, "G" * L]
        outer = ["T" * L, "A" * L]
        feat = extract_features(inner, outer)
        K = 5
        assert feat.shape == (L, K + 9)

    def test_missing_data_encoded_as_minus_one(self):
        inner = ["ANGT"]
        outer = ["NNCN"]
        feat = extract_features(inner, outer)
        # First K columns are outgroup bases; N (code 4) maps to -1.0
        assert feat[1, 0] == -1.0   # inner[0] position 1 is N
        assert feat[0, 1] == -1.0   # outer[0] position 0 is N
        assert feat[1, 1] == -1.0   # outer[0] position 1 is N

    def test_data_count_column(self):
        inner = ["ACGT"]
        outer = ["NNNN"]
        feat = extract_features(inner, outer)
        K = 2
        data_count_col = K
        # Position 0: inner has A (valid), outer has N (missing) → count=1
        assert feat[0, data_count_col] == 1.0
        # Position 3: inner has T (valid), outer has N → count=1
        assert feat[3, data_count_col] == 1.0

    def test_all_missing_data_count_zero(self):
        inner = ["NNNN"]
        outer = ["NNNN"]
        feat = extract_features(inner, outer)
        K = 2
        data_count_col = K
        assert np.all(feat[:, data_count_col] == 0.0)

    def test_agreement_ratio_perfect(self):
        inner = ["AAAA", "AAAA"]
        outer = ["AAAA"]
        feat = extract_features(inner, outer)
        K = 3
        agreement_col = K + 1
        assert np.allclose(feat[:, agreement_col], 1.0)

    def test_agreement_ratio_split(self):
        inner = ["ACGT", "TGCA"]
        outer = ["ACGT"]
        feat = extract_features(inner, outer)
        K = 3
        agreement_col = K + 1
        # At position 0: A(2), T(1) → agreement = 2/3
        assert pytest.approx(feat[0, agreement_col], abs=0.01) == 2 / 3

    def test_voting_agreement_column(self):
        inner = ["AAAA"]
        outer = ["AATT"]
        feat = extract_features(inner, outer)
        K = 2
        vote_agree_col = K + 4
        assert feat[0, vote_agree_col] == 1.0  # A == A
        assert feat[1, vote_agree_col] == 1.0  # A == A
        assert feat[2, vote_agree_col] == 0.0  # A != T
        assert feat[3, vote_agree_col] == 0.0  # A != T

    def test_with_focal_seq(self):
        inner = ["ACGT"]
        outer = ["ACGT"]
        feat_no_focal = extract_features(inner, outer)
        feat_focal = extract_features(inner, outer, focal_seq="ACGT")
        assert feat_no_focal.shape == feat_focal.shape

    def test_single_position(self):
        inner = ["A"]
        outer = ["A"]
        feat = extract_features(inner, outer)
        assert feat.shape[0] == 1


# ---------------------------------------------------------------------------
# Context helpers
# ---------------------------------------------------------------------------

class TestGCContent:

    def test_all_gc(self):
        from ancify.backend import encode_sequences
        enc = encode_sequences(["CCGGCCGG"])[0]
        gc = _gc_content(enc, window=50)
        assert np.allclose(gc, 1.0)

    def test_all_at(self):
        from ancify.backend import encode_sequences
        enc = encode_sequences(["AATTAATT"])[0]
        gc = _gc_content(enc, window=50)
        assert np.allclose(gc, 0.0)

    def test_mixed(self):
        from ancify.backend import encode_sequences
        enc = encode_sequences(["ACGT"])[0]
        gc = _gc_content(enc, window=50)
        assert pytest.approx(gc[0], abs=0.01) == 0.5

    def test_missing_excluded(self):
        from ancify.backend import encode_sequences
        enc = encode_sequences(["ACNN"])[0]
        gc = _gc_content(enc, window=50)
        # Only A and C are valid; C is GC → 1/2
        assert pytest.approx(gc[0], abs=0.01) == 0.5


class TestCpGFlag:

    def test_cg_dinucleotide(self):
        from ancify.backend import encode_sequences
        enc = encode_sequences(["ACGA"])[0]
        cpg = _cpg_flag(enc)
        assert cpg[1] == 1.0  # C at position 1
        assert cpg[2] == 1.0  # G at position 2
        assert cpg[0] == 0.0  # A before the CpG
        assert cpg[3] == 0.0  # A after the CpG

    def test_no_cg(self):
        from ancify.backend import encode_sequences
        enc = encode_sequences(["AATT"])[0]
        cpg = _cpg_flag(enc)
        assert np.all(cpg == 0.0)

    def test_single_base(self):
        from ancify.backend import encode_sequences
        enc = encode_sequences(["C"])[0]
        cpg = _cpg_flag(enc)
        assert cpg[0] == 0.0


# ---------------------------------------------------------------------------
# Training data generation
# ---------------------------------------------------------------------------

class TestLabelsFromVoting:

    def test_agreement_sites(self):
        inner = ["ACGT", "ACGT"]
        outer = ["ACGT"]
        indices, labels = _labels_from_voting(inner, outer)
        assert len(indices) == 4
        np.testing.assert_array_equal(labels, [0, 1, 2, 3])

    def test_disagreement_excluded(self):
        inner = ["AAAA"]
        outer = ["TTTT"]
        indices, labels = _labels_from_voting(inner, outer)
        assert len(indices) == 0

    def test_missing_excluded(self):
        inner = ["NNNN"]
        outer = ["AAAA"]
        indices, labels = _labels_from_voting(inner, outer)
        assert len(indices) == 0


class TestLabelsFromReference:

    def test_valid_bases(self):
        indices, labels = _labels_from_reference("ACGT")
        assert len(indices) == 4
        np.testing.assert_array_equal(labels, [0, 1, 2, 3])

    def test_n_excluded(self):
        indices, labels = _labels_from_reference("ANNG")
        assert len(indices) == 2
        np.testing.assert_array_equal(labels, [0, 2])


class TestGenerateTrainingData:

    def test_self_supervised(self):
        inner = ["ACGTACGT", "ACGTACGT"]
        outer = ["ACGTACGT"]
        feat, labels = generate_training_data(inner, outer)
        assert feat.shape[0] == labels.shape[0]
        assert feat.shape[0] == 8

    def test_reference_supervised(self):
        inner = ["ACGTACGT", "ACGTACGT"]
        outer = ["ACGTACGT"]
        feat, labels = generate_training_data(
            inner, outer, ref_seq="TGCATGCA"
        )
        assert feat.shape[0] == labels.shape[0]
        assert feat.shape[0] == 8

    def test_max_sites_caps(self):
        L = 100
        inner = ["A" * L, "A" * L]
        outer = ["A" * L]
        feat, labels = generate_training_data(inner, outer, max_sites=10)
        assert feat.shape[0] <= 10

    def test_empty_when_all_missing(self):
        inner = ["NNNN"]
        outer = ["NNNN"]
        feat, labels = generate_training_data(inner, outer)
        assert feat.shape[0] == 0


class TestStratifiedSubsample:

    def test_balances_classes(self):
        labels = np.array([0] * 100 + [1] * 100 + [2] * 100 + [3] * 100)
        rng = np.random.RandomState(42)
        keep = _stratified_subsample(labels, max_sites=40, rng=rng)
        selected_labels = labels[keep]
        for c in range(4):
            assert (selected_labels == c).sum() == 10

    def test_small_class_preserved(self):
        labels = np.array([0] * 100 + [1] * 5)
        rng = np.random.RandomState(42)
        keep = _stratified_subsample(labels, max_sites=20, rng=rng)
        selected_labels = labels[keep]
        assert (selected_labels == 1).sum() == 5


# ---------------------------------------------------------------------------
# Prediction & confidence encoding
# ---------------------------------------------------------------------------

class TestEncodeConfidence:

    def test_high_confidence_uppercase(self):
        predicted = np.array([0, 1, 2, 3], dtype=np.uint8)
        probs = np.array([0.9, 0.95, 0.85, 0.99])
        result = encode_confidence(predicted, probs)
        assert result == "ACGT"

    def test_low_confidence_lowercase(self):
        predicted = np.array([0, 1, 2, 3], dtype=np.uint8)
        probs = np.array([0.6, 0.55, 0.7, 0.5])
        result = encode_confidence(predicted, probs)
        assert result == "acgt"

    def test_unresolved_n(self):
        predicted = np.array([0, 1], dtype=np.uint8)
        probs = np.array([0.3, 0.4])
        result = encode_confidence(predicted, probs)
        assert result == "nn"

    def test_custom_thresholds(self):
        predicted = np.array([0], dtype=np.uint8)
        probs = np.array([0.6])
        assert encode_confidence(predicted, probs, high_threshold=0.5) == "A"
        assert encode_confidence(predicted, probs, high_threshold=0.9,
                                 low_threshold=0.5) == "a"
        assert encode_confidence(predicted, probs, low_threshold=0.7) == "n"

    def test_boundary_values(self):
        predicted = np.array([0, 0, 0], dtype=np.uint8)
        probs = np.array([0.8, 0.5, 0.49999])
        result = encode_confidence(predicted, probs)
        assert result[0] == "A"  # exactly at high threshold
        assert result[1] == "a"  # exactly at low threshold
        assert result[2] == "n"  # just below low threshold

    def test_empty_input(self):
        predicted = np.array([], dtype=np.uint8)
        probs = np.array([], dtype=np.float64)
        result = encode_confidence(predicted, probs)
        assert result == ""


# ---------------------------------------------------------------------------
# Model training + prediction (requires lightgbm)
# ---------------------------------------------------------------------------

FAST_PARAMS = {"num_threads": 1}


@requires_lightgbm
class TestTrainAndPredict:

    def _make_training_data(self, n=1000, rng_seed=42):
        rng = np.random.RandomState(rng_seed)
        n_features = 12
        X = rng.randn(n, n_features).astype(np.float32)
        y = rng.randint(0, 4, size=n).astype(np.int32)
        return X, y

    def test_train_returns_booster(self):
        from ancify.ml import train_model
        X, y = self._make_training_data()
        model = train_model(X, y, params=FAST_PARAMS, num_boost_round=10)
        import lightgbm as lgb
        assert isinstance(model, lgb.Booster)

    def test_predict_shape(self):
        from ancify.ml import train_model
        X, y = self._make_training_data()
        model = train_model(X, y, params=FAST_PARAMS, num_boost_round=10)
        pred, probs = predict(model, X[:10])
        assert pred.shape == (10,)
        assert probs.shape == (10,)
        assert pred.dtype == np.uint8

    def test_predict_values_valid(self):
        from ancify.ml import train_model
        X, y = self._make_training_data()
        model = train_model(X, y, params=FAST_PARAMS, num_boost_round=10)
        pred, probs = predict(model, X)
        assert np.all(pred <= 3)
        assert np.all(probs >= 0.0)
        assert np.all(probs <= 1.0)

    def test_custom_params(self):
        from ancify.ml import train_model
        X, y = self._make_training_data()
        model = train_model(
            X, y,
            params={**FAST_PARAMS, "num_leaves": 15, "learning_rate": 0.05},
            num_boost_round=10,
        )
        pred, probs = predict(model, X[:5])
        assert pred.shape == (5,)


# ---------------------------------------------------------------------------
# Model save / load round-trip (requires lightgbm)
# ---------------------------------------------------------------------------

@requires_lightgbm
class TestModelIO:

    def test_save_load_roundtrip(self, tmp_path):
        from ancify.ml import load_model, save_model, train_model
        rng = np.random.RandomState(123)
        X = rng.randn(200, 10).astype(np.float32)
        y = rng.randint(0, 4, size=200).astype(np.int32)
        model = train_model(X, y, params=FAST_PARAMS, num_boost_round=10)

        model_path = tmp_path / "test_model.lgb"
        save_model(model, str(model_path))
        assert model_path.exists()

        loaded = load_model(str(model_path))
        pred_orig, _ = predict(model, X[:20])
        pred_loaded, _ = predict(loaded, X[:20])
        np.testing.assert_array_equal(pred_orig, pred_loaded)

    def test_save_creates_parent_dirs(self, tmp_path):
        from ancify.ml import save_model, train_model
        rng = np.random.RandomState(0)
        X = rng.randn(100, 5).astype(np.float32)
        y = rng.randint(0, 4, size=100).astype(np.int32)
        model = train_model(X, y, params=FAST_PARAMS, num_boost_round=10)

        nested = tmp_path / "a" / "b" / "model.lgb"
        save_model(model, str(nested))
        assert nested.exists()


# ---------------------------------------------------------------------------
# End-to-end: extract → train → predict → encode
# ---------------------------------------------------------------------------

@requires_lightgbm
class TestEndToEnd:

    def test_pipeline_on_synthetic_data(self):
        from ancify.ml import train_model

        L = 500
        inner = ["A" * L, "A" * L, "A" * L]
        outer = ["A" * L]
        feat, labels = generate_training_data(inner, outer)

        assert np.all(labels == 0)

        rng = np.random.RandomState(99)
        extra_feat = rng.randn(100, feat.shape[1]).astype(np.float32)
        extra_labels = np.array([0] * 25 + [1] * 25 + [2] * 25 + [3] * 25,
                                dtype=np.int32)
        all_feat = np.vstack([feat, extra_feat])
        all_labels = np.concatenate([labels, extra_labels])

        model = train_model(all_feat, all_labels,
                            params=FAST_PARAMS, num_boost_round=10)
        pred, probs = predict(model, feat)
        result = encode_confidence(pred, probs)
        assert len(result) == len(feat)
        assert all(c in "ACGTacgtn" for c in result)
