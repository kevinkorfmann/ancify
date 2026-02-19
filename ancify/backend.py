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
