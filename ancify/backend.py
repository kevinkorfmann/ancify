"""Backend abstraction for CPU (NumPy) and GPU (PyTorch) compute paths.

Provides fast I/O (isal-accelerated gzip when available) and vectorized
array operations for the projection and ancestral-calling phases.
GPU paths require PyTorch with CUDA.
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)

_ORD_A = ord("A")
_ORD_C = ord("C")
_ORD_G = ord("G")
_ORD_T = ord("T")
_ORD_N = ord("N")

_UPPERCASE_LUT = np.array(
    [ord("A"), ord("C"), ord("G"), ord("T"), ord("N")], dtype=np.uint8,
)
_LOWERCASE_LUT = np.array(
    [ord("a"), ord("c"), ord("g"), ord("t"), ord("N")], dtype=np.uint8,
)

try:
    from isal import igzip
    _gz_open = igzip.open
except ImportError:
    import gzip
    _gz_open = gzip.open


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def open_gz(path):
    """Open *path* for text reading, using isal for ``.gz`` if available."""
    if str(path).endswith(".gz"):
        return _gz_open(str(path), "rt")
    return open(str(path), "rt")


# ---------------------------------------------------------------------------
# Backend / device detection
# ---------------------------------------------------------------------------

def detect_backend():
    """Return ``'gpu'`` if PyTorch with CUDA is available, else ``'cpu'``."""
    try:
        import torch
        if torch.cuda.is_available():
            return "gpu"
    except ImportError:
        pass
    return "cpu"


def get_available_gpus():
    """Return a list of usable CUDA device indices (empty if none)."""
    try:
        import torch
        if torch.cuda.is_available():
            return list(range(torch.cuda.device_count()))
    except ImportError:
        pass
    return []


def resolve_device_id(backend="auto", gpu_devices=None):
    """Choose a single GPU device ID based on configuration.

    Returns an ``int`` device index when a GPU should be used, or ``None``
    for CPU-only mode.
    """
    if backend == "cpu":
        return None
    if backend == "auto":
        backend = detect_backend()
    if backend != "gpu":
        return None
    gpus = get_available_gpus()
    if not gpus:
        return None
    if gpu_devices is not None:
        gpus = [g for g in gpus if g in gpu_devices]
    return gpus[0] if gpus else None


# ---------------------------------------------------------------------------
# Phase 1: vectorized block scatter
# ---------------------------------------------------------------------------

def _scatter_blocks_numpy(seq_array, blocks, chrom_length):
    """Fill *seq_array* from alignment *blocks* using NumPy vectorisation."""
    for start_pos, hseq_bytes, cseq_bytes in blocks:
        h = np.frombuffer(hseq_bytes, dtype=np.uint8)
        c = np.frombuffer(cseq_bytes, dtype=np.uint8)

        upper = h & 0xDF
        valid = (
            (upper == _ORD_A) | (upper == _ORD_C)
            | (upper == _ORD_G) | (upper == _ORD_T) | (upper == _ORD_N)
        )

        target = (start_pos - 1) + np.cumsum(valid, dtype=np.int64) - 1
        mask = valid & (target >= 0) & (target < chrom_length)
        seq_array[target[mask]] = c[mask]


def _scatter_blocks_gpu(seq_array, blocks, chrom_length, device_id):
    """Fill *seq_array* from alignment *blocks* on a CUDA device."""
    import torch

    device = torch.device(f"cuda:{device_id}")

    h_parts, c_parts, meta = [], [], []
    offset = 0
    for start_pos, hseq_bytes, cseq_bytes in blocks:
        n = len(hseq_bytes)
        h_parts.append(np.frombuffer(hseq_bytes, dtype=np.uint8))
        c_parts.append(np.frombuffer(cseq_bytes, dtype=np.uint8))
        meta.append((offset, start_pos, n))
        offset += n

    h_cat = torch.from_numpy(np.concatenate(h_parts)).to(device)
    c_cat = torch.from_numpy(np.concatenate(c_parts)).to(device)

    upper = h_cat & 0xDF
    valid = (
        (upper == _ORD_A) | (upper == _ORD_C)
        | (upper == _ORD_G) | (upper == _ORD_T) | (upper == _ORD_N)
    )

    target = torch.empty(offset, dtype=torch.int64, device=device)
    for cat_off, start_pos, blen in meta:
        s = slice(cat_off, cat_off + blen)
        target[s] = (
            (start_pos - 1)
            + torch.cumsum(valid[s].to(torch.int64), dim=0)
            - 1
        )

    mask = valid & (target >= 0) & (target < chrom_length)

    result = torch.from_numpy(seq_array.copy()).to(device)
    result[target[mask]] = c_cat[mask]
    np.copyto(seq_array, result.cpu().numpy())


def vectorized_block_scatter(seq_array, blocks, chrom_length, device_id=None):
    """Scatter outgroup bases from AXT *blocks* into *seq_array* (in-place).

    Parameters
    ----------
    seq_array : numpy.ndarray (uint8)
        Pre-filled array (typically all ``N``).  Modified in place.
    blocks : list of (int, bytes, bytes)
        ``(start_pos_1based, focal_seq_bytes, outgroup_seq_bytes)`` per block.
    chrom_length : int
        Length of the target chromosome (== ``len(seq_array)``).
    device_id : int or None
        CUDA device index for GPU path; ``None`` selects CPU/NumPy.
    """
    if not blocks:
        return
    if device_id is not None:
        _scatter_blocks_gpu(seq_array, blocks, chrom_length, device_id)
    else:
        _scatter_blocks_numpy(seq_array, blocks, chrom_length)


# ---------------------------------------------------------------------------
# Sequence encoding (shared by Phase 2)
# ---------------------------------------------------------------------------

def encode_sequences(seqs):
    """Encode nucleotide strings as a uint8 matrix.

    Mapping: ``A`` -> 0, ``C`` -> 1, ``G`` -> 2, ``T`` -> 3, else -> 4.
    Case-insensitive.

    Parameters
    ----------
    seqs : list of str
        Sequences of equal length.

    Returns
    -------
    numpy.ndarray, shape (n_seqs, length), dtype uint8
    """
    n = len(seqs)
    L = len(seqs[0])
    result = np.full((n, L), 4, dtype=np.uint8)
    for i, s in enumerate(seqs):
        raw = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
        upper = raw & 0xDF
        result[i, upper == _ORD_A] = 0
        result[i, upper == _ORD_C] = 1
        result[i, upper == _ORD_G] = 2
        result[i, upper == _ORD_T] = 3
    return result


# ---------------------------------------------------------------------------
# Phase 2: vectorized majority vote
# ---------------------------------------------------------------------------

def _majority_vote_numpy(encoded, min_freq):
    """CPU majority vote over columns of *encoded* (uint8 matrix)."""
    _n, L = encoded.shape
    counts = np.zeros((4, L), dtype=np.int32)
    for b in range(4):
        counts[b] = (encoded == b).sum(axis=0)
    max_count = counts.max(axis=0)
    winner = counts.argmax(axis=0).astype(np.uint8)
    winner[max_count < min_freq] = 4
    return winner


def _majority_vote_gpu_tensor(encoded_np, min_freq, device):
    """GPU majority vote; returns a ``torch.Tensor`` on *device*."""
    import torch

    enc = torch.from_numpy(encoded_np).to(device)
    _n, L = enc.shape
    counts = torch.zeros((4, L), dtype=torch.int32, device=device)
    for b in range(4):
        counts[b] = (enc == b).sum(dim=0)
    max_count, winner = counts.max(dim=0)
    winner[max_count < min_freq] = 4
    return winner.to(torch.uint8)


def vectorized_majority_vote(encoded, min_freq, device=None):
    """Vectorized column-wise majority vote.

    Parameters
    ----------
    encoded : numpy.ndarray, shape (n_seqs, L), dtype uint8
        Integer-encoded sequences (0-3 = ACGT, 4 = missing).
    min_freq : int
        Minimum count to accept a winner; positions below get 4.
    device : torch.device or None
        CUDA device for the GPU path; ``None`` uses NumPy.

    Returns
    -------
    numpy.ndarray, shape (L,), dtype uint8
        Winning base index per column (0-3), or 4 if no winner.
        Ties broken alphabetically (A > C > G > T).
    """
    if device is not None:
        return _majority_vote_gpu_tensor(encoded, min_freq, device).cpu().numpy()
    return _majority_vote_numpy(encoded, min_freq)


# ---------------------------------------------------------------------------
# Phase 2: vectorized ancestral calling
# ---------------------------------------------------------------------------

def _apply_confidence_numpy(inner_cons, outer_cons):
    """Apply confidence encoding to consensus arrays (CPU path)."""
    L = len(inner_cons)
    result = np.full(L, ord("N"), dtype=np.uint8)

    inner_valid = inner_cons < 4
    outer_valid = outer_cons < 4
    agree = inner_cons == outer_cons

    mask = inner_valid & outer_valid & agree
    result[mask] = _UPPERCASE_LUT[inner_cons[mask]]

    mask = inner_valid & ~outer_valid
    result[mask] = _LOWERCASE_LUT[inner_cons[mask]]

    mask = ~inner_valid & outer_valid
    result[mask] = _LOWERCASE_LUT[outer_cons[mask]]

    mask = inner_valid & outer_valid & ~agree
    result[mask] = ord("n")

    return result.tobytes().decode("ascii")


def _ancestral_call_gpu(inner_enc, outer_enc, min_inner, min_outer, device):
    """Full GPU pipeline: majority vote + confidence in one pass."""
    import torch

    inner_cons = _majority_vote_gpu_tensor(inner_enc, min_inner, device)
    outer_cons = _majority_vote_gpu_tensor(outer_enc, min_outer, device)

    upper_lut = torch.tensor(
        [ord("A"), ord("C"), ord("G"), ord("T"), ord("N")],
        dtype=torch.uint8, device=device,
    )
    lower_lut = torch.tensor(
        [ord("a"), ord("c"), ord("g"), ord("t"), ord("N")],
        dtype=torch.uint8, device=device,
    )

    L = inner_cons.shape[0]
    result = torch.full((L,), ord("N"), dtype=torch.uint8, device=device)

    inner_valid = inner_cons < 4
    outer_valid = outer_cons < 4
    agree = inner_cons == outer_cons

    mask = inner_valid & outer_valid & agree
    result[mask] = upper_lut[inner_cons[mask].long()]

    mask = inner_valid & ~outer_valid
    result[mask] = lower_lut[inner_cons[mask].long()]

    mask = ~inner_valid & outer_valid
    result[mask] = lower_lut[outer_cons[mask].long()]

    mask = inner_valid & outer_valid & ~agree
    result[mask] = ord("n")

    return result.cpu().numpy().tobytes().decode("ascii")


# ---------------------------------------------------------------------------
# Tree linearisation helper (shared by parsimony and likelihood GPU paths)
# ---------------------------------------------------------------------------

def _linearize_tree(tree):
    """Flatten a TreeNode into a post-order schedule for vectorised ops.

    Nodes are numbered in DFS post-order (leaves and internal nodes
    interleaved).

    Returns
    -------
    leaf_names : list of str
    leaf_indices : list of int
        Buffer index for each leaf (same order as *leaf_names*).
    ops : list of (int, list[int])
        ``(node_idx, [child_indices])`` in post-order for internal nodes.
    branch_lengths : dict[int, float]
        Node index → branch length (0.0 if absent).
    root_idx : int
    """
    leaf_names = []
    leaf_indices = []
    idx_map = {}
    ops = []
    branch_lengths = {}
    counter = [0]

    def _walk(node):
        if node.is_leaf:
            idx = counter[0]
            counter[0] += 1
            idx_map[id(node)] = idx
            leaf_names.append(node.name)
            leaf_indices.append(idx)
            branch_lengths[idx] = node.branch_length or 0.0
        else:
            for child in node.children:
                _walk(child)
            idx = counter[0]
            counter[0] += 1
            idx_map[id(node)] = idx
            child_indices = [idx_map[id(c)] for c in node.children]
            ops.append((idx, child_indices))
            branch_lengths[idx] = node.branch_length or 0.0

    _walk(tree)
    root_idx = idx_map[id(tree)]
    return leaf_names, leaf_indices, ops, branch_lengths, root_idx


# ---------------------------------------------------------------------------
# Bitmask encoding for Fitch parsimony: A=1, C=2, G=4, T=8, N=15
# ---------------------------------------------------------------------------

_BITMASK_LUT_NP = np.array([1, 2, 4, 8, 15], dtype=np.uint8)

_LOWEST_BIT_LUT_NP = np.array(
    [4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0], dtype=np.uint8,
)

_POPCOUNT_LUT_NP = np.array(
    [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4], dtype=np.uint8,
)


# ---------------------------------------------------------------------------
# Phase 2: vectorised Fitch parsimony
# ---------------------------------------------------------------------------

def _fitch_parsimony_numpy(species_seqs, tree):
    """CPU (NumPy) vectorised Fitch parsimony for an entire chromosome."""
    leaf_names, leaf_indices, ops, _, root_idx = _linearize_tree(tree)
    L = len(species_seqs[leaf_names[0]])
    n_nodes = root_idx + 1

    encoded = encode_sequences([species_seqs[n] for n in leaf_names])
    leaf_bitmasks = _BITMASK_LUT_NP[encoded]

    buf = np.zeros((n_nodes, L), dtype=np.uint8)
    for i, li in enumerate(leaf_indices):
        buf[li] = leaf_bitmasks[i]

    for node_idx, children in ops:
        intersection = buf[children[0]]
        for ci in children[1:]:
            intersection = intersection & buf[ci]
        union = buf[children[0]]
        for ci in children[1:]:
            union = union | buf[ci]
        buf[node_idx] = np.where(intersection > 0, intersection, union)

    root_bm = buf[root_idx]
    root_allele = _LOWEST_BIT_LUT_NP[root_bm]
    popcount = _POPCOUNT_LUT_NP[root_bm]

    any_data = np.zeros(L, dtype=bool)
    for i in range(len(leaf_names)):
        any_data |= leaf_bitmasks[i] != 15

    result = np.full(L, ord("N"), dtype=np.uint8)
    mask = popcount == 1
    result[mask] = _UPPERCASE_LUT[root_allele[mask]]
    mask = (popcount > 1) & (popcount < 4)
    result[mask] = _LOWERCASE_LUT[root_allele[mask]]
    mask = (root_bm == 15) & any_data
    result[mask] = ord("n")

    return result.tobytes().decode("ascii")


def _fitch_parsimony_gpu(species_seqs, tree, device):
    """GPU (PyTorch) Fitch parsimony for an entire chromosome."""
    import torch

    leaf_names, leaf_indices, ops, _, root_idx = _linearize_tree(tree)
    L = len(species_seqs[leaf_names[0]])
    n_nodes = root_idx + 1

    encoded = encode_sequences([species_seqs[n] for n in leaf_names])
    bitmask_lut = torch.tensor([1, 2, 4, 8, 15], dtype=torch.uint8, device=device)
    leaf_bitmasks = bitmask_lut[torch.from_numpy(encoded).long().to(device)]

    buf = torch.zeros((n_nodes, L), dtype=torch.uint8, device=device)
    for i, li in enumerate(leaf_indices):
        buf[li] = leaf_bitmasks[i]

    for node_idx, children in ops:
        intersection = buf[children[0]]
        for ci in children[1:]:
            intersection = intersection & buf[ci]
        union_val = buf[children[0]]
        for ci in children[1:]:
            union_val = union_val | buf[ci]
        buf[node_idx] = torch.where(intersection > 0, intersection, union_val)

    root_bm = buf[root_idx]

    lowest_bit_lut = torch.tensor(
        [4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0],
        dtype=torch.uint8, device=device,
    )
    popcount_lut = torch.tensor(
        [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4],
        dtype=torch.uint8, device=device,
    )
    root_allele = lowest_bit_lut[root_bm.long()]
    popcount = popcount_lut[root_bm.long()]

    any_data = torch.zeros(L, dtype=torch.bool, device=device)
    for i in range(len(leaf_names)):
        any_data |= leaf_bitmasks[i] != 15

    upper_lut = torch.tensor(
        [ord("A"), ord("C"), ord("G"), ord("T"), ord("N")],
        dtype=torch.uint8, device=device,
    )
    lower_lut = torch.tensor(
        [ord("a"), ord("c"), ord("g"), ord("t"), ord("N")],
        dtype=torch.uint8, device=device,
    )

    result = torch.full((L,), ord("N"), dtype=torch.uint8, device=device)
    mask = popcount == 1
    result[mask] = upper_lut[root_allele[mask].long()]
    mask = (popcount > 1) & (popcount < 4)
    result[mask] = lower_lut[root_allele[mask].long()]
    mask = (root_bm == 15) & any_data
    result[mask] = ord("n")

    return result.cpu().numpy().tobytes().decode("ascii")


def vectorized_fitch_call(species_seqs, tree, device=None):
    """Vectorised Fitch parsimony ancestral calling for a chromosome.

    Parameters
    ----------
    species_seqs : dict[str, str]
        Mapping of species name to projected sequence (same length).
    tree : TreeNode
        Phylogenetic tree whose leaves match keys in *species_seqs*.
    device : torch.device or None
        CUDA device for the GPU path; ``None`` uses NumPy.

    Returns
    -------
    str
        Ancestral sequence with case-encoded confidence.
    """
    if device is not None:
        return _fitch_parsimony_gpu(species_seqs, tree, device)
    return _fitch_parsimony_numpy(species_seqs, tree)


# ---------------------------------------------------------------------------
# Phase 2: vectorised Felsenstein likelihood
# ---------------------------------------------------------------------------

def _leaf_likelihoods(encoded_row, device=None):
    """Convert an encoded leaf row to a (4, L) likelihood array/tensor.

    Known bases get one-hot columns; missing (4) gets all-ones.
    """
    L = encoded_row.shape[0]
    if device is not None:
        import torch
        enc = torch.from_numpy(encoded_row).long().to(device)
        lik = torch.zeros((4, L), dtype=torch.float32, device=device)
        for b in range(4):
            lik[b, enc == b] = 1.0
        lik[:, enc == 4] = 1.0
        return lik
    else:
        enc = encoded_row
        lik = np.zeros((4, L), dtype=np.float32)
        for b in range(4):
            lik[b, enc == b] = 1.0
        lik[:, enc == 4] = 1.0
        return lik


def _felsenstein_numpy(species_seqs, tree, model, high_thresh, low_thresh):
    """CPU (NumPy) vectorised Felsenstein pruning for an entire chromosome."""
    leaf_names, leaf_indices, ops, branch_lengths, root_idx = _linearize_tree(tree)
    L = len(species_seqs[leaf_names[0]])
    n_nodes = root_idx + 1

    encoded = encode_sequences([species_seqs[n] for n in leaf_names])

    P_mats = {}
    for idx, bl in branch_lengths.items():
        P_mats[idx] = model.transition_probs(bl)

    lik_buf = [None] * n_nodes
    for i, li in enumerate(leaf_indices):
        lik_buf[li] = _leaf_likelihoods(encoded[i])

    for node_idx, children in ops:
        node_lik = np.ones((4, L), dtype=np.float32)
        for ci in children:
            P = P_mats[ci]
            node_lik *= P @ lik_buf[ci]
        lik_buf[node_idx] = node_lik

    pi = model.base_frequencies()
    joint = pi[:, None] * lik_buf[root_idx]
    total = joint.sum(axis=0, keepdims=True)
    total = np.where(total == 0, 1.0, total)
    posterior = joint / total

    max_prob = posterior.max(axis=0)
    best_base = posterior.argmax(axis=0).astype(np.uint8)

    any_data = np.zeros(L, dtype=bool)
    for i in range(len(leaf_names)):
        any_data |= encoded[i] != 4

    result = np.full(L, ord("N"), dtype=np.uint8)
    mask = any_data & (max_prob >= high_thresh)
    result[mask] = _UPPERCASE_LUT[best_base[mask]]
    mask = any_data & (max_prob >= low_thresh) & (max_prob < high_thresh)
    result[mask] = _LOWERCASE_LUT[best_base[mask]]
    mask = any_data & (max_prob < low_thresh)
    result[mask] = ord("n")

    return result.tobytes().decode("ascii")


def _felsenstein_gpu(species_seqs, tree, model, high_thresh, low_thresh, device):
    """GPU (PyTorch) Felsenstein pruning for an entire chromosome."""
    import torch

    leaf_names, leaf_indices, ops, branch_lengths, root_idx = _linearize_tree(tree)
    L = len(species_seqs[leaf_names[0]])
    n_nodes = root_idx + 1

    encoded = encode_sequences([species_seqs[n] for n in leaf_names])

    P_tensors = {}
    for idx, bl in branch_lengths.items():
        P_np = model.transition_probs(bl)
        P_tensors[idx] = torch.from_numpy(P_np.astype(np.float32)).to(device)

    lik_buf = [None] * n_nodes
    for i, li in enumerate(leaf_indices):
        lik_buf[li] = _leaf_likelihoods(encoded[i], device=device)

    for node_idx, children in ops:
        node_lik = torch.ones((4, L), dtype=torch.float32, device=device)
        for ci in children:
            P = P_tensors[ci]
            node_lik = node_lik * (P @ lik_buf[ci])
        lik_buf[node_idx] = node_lik

    pi = torch.from_numpy(
        model.base_frequencies().astype(np.float32),
    ).to(device)
    joint = pi[:, None] * lik_buf[root_idx]
    total = joint.sum(dim=0, keepdim=True)
    total = torch.where(total == 0, torch.ones_like(total), total)
    posterior = joint / total

    max_prob, best_base = posterior.max(dim=0)
    best_base = best_base.to(torch.uint8)

    enc_t = torch.from_numpy(encoded).to(device)
    any_data = (enc_t != 4).any(dim=0)

    upper_lut = torch.tensor(
        [ord("A"), ord("C"), ord("G"), ord("T"), ord("N")],
        dtype=torch.uint8, device=device,
    )
    lower_lut = torch.tensor(
        [ord("a"), ord("c"), ord("g"), ord("t"), ord("N")],
        dtype=torch.uint8, device=device,
    )

    result = torch.full((L,), ord("N"), dtype=torch.uint8, device=device)
    mask = any_data & (max_prob >= high_thresh)
    result[mask] = upper_lut[best_base[mask].long()]
    mask = any_data & (max_prob >= low_thresh) & (max_prob < high_thresh)
    result[mask] = lower_lut[best_base[mask].long()]
    mask = any_data & (max_prob < low_thresh)
    result[mask] = ord("n")

    return result.cpu().numpy().tobytes().decode("ascii")


def vectorized_likelihood_call(
    species_seqs, tree, model,
    high_threshold=0.8, low_threshold=0.5,
    device=None,
):
    """Vectorised Felsenstein likelihood ancestral calling for a chromosome.

    Parameters
    ----------
    species_seqs : dict[str, str]
        Mapping of species name to projected sequence (same length).
    tree : TreeNode
        Phylogenetic tree whose leaves match keys in *species_seqs*.
    model : SubstitutionModel
        Instantiated substitution model.
    high_threshold, low_threshold : float
        Posterior probability thresholds for confidence encoding.
    device : torch.device or None
        CUDA device for the GPU path; ``None`` uses NumPy.

    Returns
    -------
    str
        Ancestral sequence with case-encoded confidence.
    """
    if device is not None:
        return _felsenstein_gpu(
            species_seqs, tree, model,
            high_threshold, low_threshold, device,
        )
    return _felsenstein_numpy(
        species_seqs, tree, model,
        high_threshold, low_threshold,
    )


# ---------------------------------------------------------------------------
# Phase 2: vectorized ancestral calling (voting)
# ---------------------------------------------------------------------------

def vectorized_ancestral_call(
    inner_seqs, outer_seqs,
    min_inner_freq=1, min_outer_freq=1,
    device=None,
):
    """Vectorized ancestral allele calling for an entire chromosome.

    Encodes sequences, runs majority vote on each outgroup tier, then
    assigns confidence-coded bases in a single vectorized pass.

    Parameters
    ----------
    inner_seqs, outer_seqs : list of str
        Projected outgroup sequences (same length).
    min_inner_freq, min_outer_freq : int
        Minimum allele count for the majority vote in each tier.
    device : torch.device or None
        CUDA device for the GPU path; ``None`` uses NumPy.

    Returns
    -------
    str
        Ancestral sequence with case-encoded confidence.
    """
    inner_enc = encode_sequences(inner_seqs)
    outer_enc = encode_sequences(outer_seqs)

    if device is not None:
        return _ancestral_call_gpu(
            inner_enc, outer_enc, min_inner_freq, min_outer_freq, device,
        )

    inner_cons = _majority_vote_numpy(inner_enc, min_inner_freq)
    outer_cons = _majority_vote_numpy(outer_enc, min_outer_freq)
    return _apply_confidence_numpy(inner_cons, outer_cons)
