"""Machine learning-based ancestral allele calling.

Uses a LightGBM gradient-boosted classifier trained on per-position
features extracted from outgroup alignments.  Confidence is calibrated
from predicted class probabilities rather than the binary agree/disagree
scheme of the voting method.

Optional dependency: requires ``lightgbm`` and ``scikit-learn``.
Install with ``pip install ancify[ml]``.
"""

import logging
from pathlib import Path

import numpy as np

from .backend import encode_sequences

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Context helpers
# ---------------------------------------------------------------------------

def _gc_content(seq_encoded, window):
    """GC fraction in a sliding window of half-width *window*."""
    L = len(seq_encoded)
    gc = ((seq_encoded == 1) | (seq_encoded == 2)).astype(np.float64)
    valid = (seq_encoded < 4).astype(np.float64)
    w = min(window, (L - 1) // 2) if L > 1 else 0
    if w == 0:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(valid > 0, gc / valid, 0.0).astype(np.float32)
    kernel = np.ones(2 * w + 1, dtype=np.float64)
    gc_sum = np.convolve(gc, kernel, mode="same")
    valid_sum = np.convolve(valid, kernel, mode="same")
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(valid_sum > 0, gc_sum / valid_sum, 0.0).astype(np.float32)


def _cpg_flag(consensus_encoded):
    """Boolean flag: position is part of a CpG dinucleotide."""
    L = len(consensus_encoded)
    cpg = np.zeros(L, dtype=bool)
    if L > 1:
        cg = (consensus_encoded[:-1] == 1) & (consensus_encoded[1:] == 2)
        cpg[:-1] |= cg
        cpg[1:] |= cg
    return cpg.astype(np.float32)


# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------

def extract_features(inner_seqs, outer_seqs, focal_seq=None, window=50):
    """Build a per-position feature matrix from projected outgroup sequences.

    Parameters
    ----------
    inner_seqs : list of str
        Projected inner outgroup sequences (equal length).
    outer_seqs : list of str
        Projected outer outgroup sequences (equal length).
    focal_seq : str or None
        Focal reference sequence for context features.  When *None*,
        context is derived from the outgroup majority consensus.
    window : int
        Half-window size for GC content (default 50 bp each side).

    Returns
    -------
    numpy.ndarray, shape ``(L, n_features)``, dtype float32
    """
    all_seqs = inner_seqs + outer_seqs
    encoded = encode_sequences(all_seqs)          # (K, L) uint8
    K, L = encoded.shape
    n_inner = len(inner_seqs)

    features = []

    # 1) Outgroup bases: 0-3 for ACGT, -1 for missing
    outgroup_feat = encoded.astype(np.float32)
    outgroup_feat[encoded == 4] = -1.0
    features.append(outgroup_feat.T)              # (L, K)

    # 2) Data count
    valid_mask = encoded < 4                      # (K, L)
    data_count = valid_mask.sum(axis=0).astype(np.float32)
    features.append(data_count[:, np.newaxis])    # (L, 1)

    # 3) Agreement ratio
    counts = np.zeros((4, L), dtype=np.int32)
    for b in range(4):
        counts[b] = (encoded == b).sum(axis=0)
    max_count = counts.max(axis=0).astype(np.float32)
    valid_total = data_count
    with np.errstate(divide="ignore", invalid="ignore"):
        agreement = np.where(valid_total > 0, max_count / valid_total, 0.0)
    features.append(agreement[:, np.newaxis])     # (L, 1)

    # 4) Inner / outer consensus (0-4)
    inner_enc = encoded[:n_inner]
    inner_counts = np.zeros((4, L), dtype=np.int32)
    for b in range(4):
        inner_counts[b] = (inner_enc == b).sum(axis=0)
    inner_cons = inner_counts.argmax(axis=0).astype(np.float32)
    inner_cons[inner_counts.max(axis=0) == 0] = 4.0

    outer_enc = encoded[n_inner:]
    outer_counts = np.zeros((4, L), dtype=np.int32)
    for b in range(4):
        outer_counts[b] = (outer_enc == b).sum(axis=0)
    outer_cons = outer_counts.argmax(axis=0).astype(np.float32)
    outer_cons[outer_counts.max(axis=0) == 0] = 4.0

    features.append(inner_cons[:, np.newaxis])    # (L, 1)
    features.append(outer_cons[:, np.newaxis])    # (L, 1)

    # 5) Voting agreement (inner == outer, both valid)
    vote_agree = ((inner_cons == outer_cons)
                  & (inner_cons < 4)
                  & (outer_cons < 4)).astype(np.float32)
    features.append(vote_agree[:, np.newaxis])    # (L, 1)

    # 6-9) Context features from focal ref or outgroup consensus
    if focal_seq is not None:
        ctx = encode_sequences([focal_seq])[0]
    else:
        ctx = counts.argmax(axis=0).astype(np.uint8)
        ctx[counts.max(axis=0) == 0] = 4

    features.append(_gc_content(ctx, window)[:, np.newaxis])   # (L, 1)
    features.append(_cpg_flag(ctx)[:, np.newaxis])             # (L, 1)

    flank_up = np.full(L, 4.0, dtype=np.float32)
    flank_down = np.full(L, 4.0, dtype=np.float32)
    if L > 1:
        flank_up[1:] = ctx[:-1].astype(np.float32)
        flank_down[:-1] = ctx[1:].astype(np.float32)
    features.append(flank_up[:, np.newaxis])                   # (L, 1)
    features.append(flank_down[:, np.newaxis])                 # (L, 1)

    return np.hstack(features).astype(np.float32)


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------

def train_model(features, labels, params=None, num_boost_round=200):
    """Train a LightGBM multi-class classifier.

    Parameters
    ----------
    features : numpy.ndarray, shape ``(N, n_features)``
    labels : numpy.ndarray, shape ``(N,)``, values in 0-3 (A/C/G/T)
    params : dict or None
        LightGBM training parameters; merged over sensible defaults.
    num_boost_round : int
        Number of boosting iterations.

    Returns
    -------
    lightgbm.Booster
    """
    import lightgbm as lgb

    defaults = {
        "objective": "multiclass",
        "num_class": 4,
        "metric": "multi_logloss",
        "learning_rate": 0.1,
        "num_leaves": 63,
        "min_data_in_leaf": 100,
        "feature_fraction": 0.8,
        "bagging_fraction": 0.8,
        "bagging_freq": 5,
        "verbose": -1,
    }
    if params is not None:
        defaults.update(params)

    dataset = lgb.Dataset(features, label=labels)
    return lgb.train(defaults, dataset, num_boost_round=num_boost_round)


# ---------------------------------------------------------------------------
# Prediction & confidence encoding
# ---------------------------------------------------------------------------

def predict(model, features):
    """Run inference on a feature matrix.

    Returns
    -------
    predicted : numpy.ndarray, shape ``(N,)``, dtype uint8
        Predicted base index (0-3).
    probabilities : numpy.ndarray, shape ``(N,)``
        Probability of the winning class.
    """
    probs = model.predict(features)               # (N, 4)
    predicted = probs.argmax(axis=1).astype(np.uint8)
    probabilities = probs.max(axis=1)
    return predicted, probabilities


def encode_confidence(predicted, probabilities,
                      high_threshold=0.8, low_threshold=0.5):
    """Map predictions + probabilities to confidence-coded characters.

    * ``p >= high_threshold`` --> uppercase (high confidence)
    * ``p >= low_threshold``  --> lowercase (low confidence)
    * ``p <  low_threshold``  --> ``n`` (unresolved)

    Positions where all outgroups are missing should be set to ``N``
    by the caller after this function returns.

    Returns
    -------
    str
        Confidence-encoded sequence of length ``len(predicted)``.
    """
    upper = np.array([ord("A"), ord("C"), ord("G"), ord("T")], dtype=np.uint8)
    lower = np.array([ord("a"), ord("c"), ord("g"), ord("t")], dtype=np.uint8)

    result = np.full(len(predicted), ord("n"), dtype=np.uint8)
    hi = probabilities >= high_threshold
    lo = (probabilities >= low_threshold) & ~hi
    result[hi] = upper[predicted[hi]]
    result[lo] = lower[predicted[lo]]
    return result.tobytes().decode("ascii")


# ---------------------------------------------------------------------------
# Model I/O
# ---------------------------------------------------------------------------

def save_model(model, path):
    """Persist a trained LightGBM Booster to *path*."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    model.save_model(str(path))
    logger.info("Model saved to %s", path)


def load_model(path):
    """Load a LightGBM Booster from *path*."""
    import lightgbm as lgb

    model = lgb.Booster(model_file=str(path))
    logger.info("Model loaded from %s", path)
    return model


# ---------------------------------------------------------------------------
# Training data generation
# ---------------------------------------------------------------------------

def _labels_from_voting(inner_seqs, outer_seqs):
    """Derive training labels from high-confidence voting sites.

    Returns arrays of indices and labels (0-3) for positions where
    inner and outer consensus agree (uppercase in voting output).
    """
    inner_enc = encode_sequences(inner_seqs)
    outer_enc = encode_sequences(outer_seqs)
    L = inner_enc.shape[1]

    inner_counts = np.zeros((4, L), dtype=np.int32)
    outer_counts = np.zeros((4, L), dtype=np.int32)
    for b in range(4):
        inner_counts[b] = (inner_enc == b).sum(axis=0)
        outer_counts[b] = (outer_enc == b).sum(axis=0)

    inner_winner = inner_counts.argmax(axis=0)
    outer_winner = outer_counts.argmax(axis=0)
    inner_valid = inner_counts.max(axis=0) > 0
    outer_valid = outer_counts.max(axis=0) > 0
    agree = inner_valid & outer_valid & (inner_winner == outer_winner)

    indices = np.where(agree)[0]
    labels = inner_winner[agree].astype(np.int32)
    return indices, labels


def _labels_from_reference(ref_seq):
    """Derive training labels from an external reference FASTA.

    Returns indices and labels for positions where the reference has
    a valid base (A/C/G/T).
    """
    ref_enc = encode_sequences([ref_seq])[0]
    valid = ref_enc < 4
    indices = np.where(valid)[0]
    labels = ref_enc[valid].astype(np.int32)
    return indices, labels


def generate_training_data(inner_seqs, outer_seqs, ref_seq=None,
                           max_sites=5_000_000, rng_seed=42):
    """Build a subsampled training dataset for one chromosome.

    Parameters
    ----------
    inner_seqs, outer_seqs : list of str
        Projected outgroup sequences.
    ref_seq : str or None
        External reference ancestral sequence.  If *None*, labels are
        derived from high-confidence voting sites (self-supervised).
    max_sites : int
        Maximum training sites to return (stratified subsample).
    rng_seed : int
        Random seed for reproducible subsampling.

    Returns
    -------
    features : numpy.ndarray
    labels : numpy.ndarray
    """
    if ref_seq is not None:
        indices, labels = _labels_from_reference(ref_seq)
    else:
        indices, labels = _labels_from_voting(inner_seqs, outer_seqs)

    if len(indices) == 0:
        return np.empty((0, 0), dtype=np.float32), np.empty(0, dtype=np.int32)

    feat_all = extract_features(inner_seqs, outer_seqs)

    if len(indices) > max_sites:
        rng = np.random.RandomState(rng_seed)
        keep = _stratified_subsample(labels, max_sites, rng)
        indices = indices[keep]
        labels = labels[keep]

    return feat_all[indices], labels


def _stratified_subsample(labels, max_sites, rng):
    """Subsample indices with equal representation per class."""
    classes = np.unique(labels)
    per_class = max_sites // len(classes)
    selected = []
    for c in classes:
        idx = np.where(labels == c)[0]
        if len(idx) > per_class:
            idx = rng.choice(idx, per_class, replace=False)
        selected.append(idx)
    result = np.concatenate(selected)
    rng.shuffle(result)
    return result


def train_from_config(config):
    """End-to-end training workflow driven by a PipelineConfig.

    Iterates over chromosomes, extracts training data, trains a single
    LightGBM model, and saves it to ``config.ml_model_path``.
    """
    from .utils import read_fasta

    chromosomes = config.resolve_chromosomes()
    work = Path(config.work_dir)

    all_features, all_labels = [], []
    ref_dir = getattr(config, "ml_training_reference", None)

    for chrom in chromosomes:
        inner_seqs = [
            read_fasta(str(work / "projected" / og.name / f"{chrom}.fa"))[1]
            for og in config.outgroups_inner
        ]
        outer_seqs = [
            read_fasta(str(work / "projected" / og.name / f"{chrom}.fa"))[1]
            for og in config.outgroups_outer
        ]

        ref_seq = None
        if ref_dir:
            ref_path = Path(ref_dir) / f"{chrom}.fa"
            if ref_path.exists():
                _, ref_seq = read_fasta(str(ref_path))

        feat, lab = generate_training_data(
            inner_seqs, outer_seqs, ref_seq=ref_seq,
        )
        if feat.size:
            all_features.append(feat)
            all_labels.append(lab)
        logger.info("  Training data from %s: %d sites", chrom, len(lab))

    if not all_features:
        raise RuntimeError("No training data could be generated from any chromosome.")

    features = np.vstack(all_features)
    labels = np.concatenate(all_labels)
    logger.info("Training on %d total sites", len(labels))

    model = train_model(features, labels)

    out_path = config.ml_model_path
    if out_path is None:
        out_path = str(Path(config.work_dir) / "ancify_ml_model.lgb")
    save_model(model, out_path)
    return out_path


# ---------------------------------------------------------------------------
# Per-chromosome worker (parallel dispatch target)
# ---------------------------------------------------------------------------

def _call_chromosome_ml(args):
    """Worker: call ancestral states for one chromosome using ML.

    Argument tuple matches the dispatch pattern of
    ``_call_chromosome_vectorized`` in :mod:`ancify.ancestral`.
    """
    (chrom, inner_paths, outer_paths, out_path,
     model_path, high_threshold, low_threshold) = args

    from .utils import read_fasta, write_fasta

    inner_seqs = [read_fasta(p)[1] for p in inner_paths]
    outer_seqs = [read_fasta(p)[1] for p in outer_paths]

    length = len(inner_seqs[0])
    for s in inner_seqs + outer_seqs:
        if len(s) != length:
            raise ValueError(
                f"Length mismatch on {chrom}: expected {length}, got {len(s)}"
            )

    encoded = encode_sequences(inner_seqs + outer_seqs)
    all_missing = (encoded == 4).all(axis=0)

    features = extract_features(inner_seqs, outer_seqs)
    model = load_model(model_path)
    predicted, probabilities = predict(model, features)

    anc = encode_confidence(predicted, probabilities, high_threshold, low_threshold)
    anc_arr = np.frombuffer(anc.encode("ascii"), dtype=np.uint8).copy()
    anc_arr[all_missing] = ord("N")

    write_fasta(out_path, f">{chrom}", anc_arr.tobytes().decode("ascii"))
    return chrom
