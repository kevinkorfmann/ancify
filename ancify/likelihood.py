"""Likelihood-based ancestral reconstruction via Felsenstein's pruning algorithm.

Implements four continuous-time Markov substitution models (JC69, K80, HKY85,
GTR) and Felsenstein's pruning algorithm to compute posterior probabilities of
ancestral bases at the root of a phylogenetic tree.

Each model defines a 4×4 instantaneous rate matrix Q normalised so that the
expected number of substitutions per unit branch length is 1.  Transition
probability matrices P(t) = expm(Q·t) are computed via ``scipy.linalg.expm``
(or a closed-form shortcut for JC69).
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from typing import Dict, List, Optional

import numpy as np
from scipy.linalg import expm

from .parsimony import TreeNode, get_leaf_names, parse_newick
from .utils import read_fasta, write_fasta

BASES = ["A", "C", "G", "T"]
BASE_INDEX = {b: i for i, b in enumerate(BASES)}

# Purine ↔ purine (A↔G) and pyrimidine ↔ pyrimidine (C↔T) index pairs.
_TRANSITIONS = frozenset({(0, 2), (2, 0), (1, 3), (3, 1)})


# ---------------------------------------------------------------------------
# Substitution models
# ---------------------------------------------------------------------------

def _normalize_Q(Q: np.ndarray, pi: np.ndarray) -> np.ndarray:
    """Scale *Q* in-place so the expected substitution rate equals 1."""
    rate = -float(np.dot(pi, np.diag(Q)))
    if rate > 0:
        Q /= rate
    return Q


class SubstitutionModel(ABC):
    """Base class for continuous-time Markov substitution models."""

    @abstractmethod
    def rate_matrix(self) -> np.ndarray:
        """Return the normalised 4×4 instantaneous rate matrix Q."""
        ...

    @abstractmethod
    def base_frequencies(self) -> np.ndarray:
        """Return equilibrium frequencies [π_A, π_C, π_G, π_T]."""
        ...

    def transition_probs(self, t: float) -> np.ndarray:
        """Return 4×4 transition probability matrix P(t) = expm(Q·t)."""
        if t <= 0.0:
            return np.eye(4)
        return expm(self.rate_matrix() * t)


class JC69(SubstitutionModel):
    """Jukes–Cantor (1969): all substitutions equally likely."""

    def base_frequencies(self) -> np.ndarray:
        return np.full(4, 0.25)

    def rate_matrix(self) -> np.ndarray:
        Q = np.ones((4, 4))
        np.fill_diagonal(Q, 0.0)
        for i in range(4):
            Q[i, i] = -Q[i].sum()
        return _normalize_Q(Q, self.base_frequencies())

    def transition_probs(self, t: float) -> np.ndarray:
        if t <= 0.0:
            return np.eye(4)
        e = math.exp(-4.0 * t / 3.0)
        p_same = 0.25 + 0.75 * e
        p_diff = 0.25 - 0.25 * e
        P = np.full((4, 4), p_diff)
        np.fill_diagonal(P, p_same)
        return P


class K80(SubstitutionModel):
    """Kimura (1980) two-parameter model with transition/transversion ratio κ."""

    def __init__(self, kappa: float = 2.0):
        self.kappa = kappa

    def base_frequencies(self) -> np.ndarray:
        return np.full(4, 0.25)

    def rate_matrix(self) -> np.ndarray:
        Q = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                if i != j:
                    Q[i, j] = self.kappa if (i, j) in _TRANSITIONS else 1.0
            Q[i, i] = -Q[i].sum()
        return _normalize_Q(Q, self.base_frequencies())


class HKY85(SubstitutionModel):
    """Hasegawa–Kishino–Yano (1985): κ with non-uniform base frequencies."""

    def __init__(self, kappa: float = 2.0, pi: Optional[List[float]] = None):
        self.kappa = kappa
        self._pi = np.asarray(pi if pi is not None else [0.25] * 4, dtype=float)

    def base_frequencies(self) -> np.ndarray:
        return self._pi.copy()

    def rate_matrix(self) -> np.ndarray:
        pi = self._pi
        Q = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                if i != j:
                    Q[i, j] = (self.kappa if (i, j) in _TRANSITIONS else 1.0) * pi[j]
            Q[i, i] = -Q[i].sum()
        return _normalize_Q(Q, pi)


class GTR(SubstitutionModel):
    """General time-reversible model: 6 exchangeability rates + base freqs.

    *rates* order: [AC, AG, AT, CG, CT, GT].
    """

    _PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

    def __init__(
        self,
        rates: Optional[List[float]] = None,
        pi: Optional[List[float]] = None,
    ):
        self._rates = list(rates) if rates is not None else [1.0] * 6
        self._pi = np.asarray(pi if pi is not None else [0.25] * 4, dtype=float)

    def base_frequencies(self) -> np.ndarray:
        return self._pi.copy()

    def rate_matrix(self) -> np.ndarray:
        pi = self._pi
        Q = np.zeros((4, 4))
        for idx, (i, j) in enumerate(self._PAIRS):
            Q[i, j] = self._rates[idx] * pi[j]
            Q[j, i] = self._rates[idx] * pi[i]
        for i in range(4):
            Q[i, i] = -Q[i].sum()
        return _normalize_Q(Q, pi)


# ---------------------------------------------------------------------------
# Model factory
# ---------------------------------------------------------------------------

def build_model(
    name: str,
    kappa: float = 2.0,
    base_freqs: Optional[List[float]] = None,
    rates: Optional[List[float]] = None,
) -> SubstitutionModel:
    """Instantiate a substitution model by name."""
    upper = name.upper()
    if upper == "JC69":
        return JC69()
    if upper == "K80":
        return K80(kappa=kappa)
    if upper == "HKY85":
        return HKY85(kappa=kappa, pi=base_freqs)
    if upper == "GTR":
        return GTR(rates=rates, pi=base_freqs)
    raise ValueError(f"Unknown substitution model: {name!r}")


# ---------------------------------------------------------------------------
# Felsenstein pruning
# ---------------------------------------------------------------------------

def _build_P_cache(tree: TreeNode, model: SubstitutionModel) -> Dict[int, np.ndarray]:
    """Pre-compute P(t) for every edge in *tree*."""
    cache: Dict[int, np.ndarray] = {}

    def _walk(node: TreeNode) -> None:
        t = node.branch_length if node.branch_length is not None else 0.0
        cache[id(node)] = model.transition_probs(t)
        for child in node.children:
            _walk(child)

    _walk(tree)
    return cache


def _pruning_up(
    node: TreeNode,
    leaf_alleles: Dict[str, str],
    P_cache: Dict[int, np.ndarray],
) -> np.ndarray:
    """Recursive bottom-up pass returning the conditional likelihood vector."""
    if node.is_leaf:
        allele = leaf_alleles.get(node.name, "N").upper()
        if allele in BASE_INDEX:
            lik = np.zeros(4)
            lik[BASE_INDEX[allele]] = 1.0
        else:
            lik = np.ones(4)
        return lik

    lik = np.ones(4)
    for child in node.children:
        child_lik = _pruning_up(child, leaf_alleles, P_cache)
        P = P_cache[id(child)]
        lik *= P @ child_lik
    return lik


def felsenstein_pruning(
    tree: TreeNode,
    leaf_alleles: Dict[str, str],
    model: SubstitutionModel,
    _P_cache: Optional[Dict[int, np.ndarray]] = None,
) -> Dict[str, float]:
    """Posterior probabilities for the root base via Felsenstein's algorithm.

    Returns
    -------
    dict
        ``{"A": p, "C": p, "G": p, "T": p}`` summing to 1.
    """
    if _P_cache is None:
        _P_cache = _build_P_cache(tree, model)

    root_lik = _pruning_up(tree, leaf_alleles, _P_cache)
    pi = model.base_frequencies()
    joint = pi * root_lik
    total = joint.sum()

    if total == 0.0:
        return {b: 0.25 for b in BASES}

    posterior = joint / total
    return {b: float(posterior[i]) for i, b in enumerate(BASES)}


# ---------------------------------------------------------------------------
# Per-position ancestral caller
# ---------------------------------------------------------------------------

def call_ancestral_base_likelihood(
    tree: TreeNode,
    species_bases: Dict[str, str],
    model: SubstitutionModel,
    high_threshold: float = 0.8,
    low_threshold: float = 0.5,
    _P_cache: Optional[Dict[int, np.ndarray]] = None,
) -> str:
    """Infer the ancestral allele at one position using likelihood.

    Returns
    -------
    str
        Uppercase (high confidence), lowercase (low confidence),
        ``'n'`` (unresolved), or ``'N'`` (all missing).
    """
    leaf_names = get_leaf_names(tree)
    has_data = any(
        species_bases.get(name, "N").upper() in BASE_INDEX
        for name in leaf_names
    )
    if not has_data:
        return "N"

    posterior = felsenstein_pruning(tree, species_bases, model, _P_cache=_P_cache)
    best_base = max(posterior, key=posterior.get)
    p = posterior[best_base]

    if p >= high_threshold:
        return best_base.upper()
    if p >= low_threshold:
        return best_base.lower()
    return "n"


# ---------------------------------------------------------------------------
# Per-chromosome parallel worker
# ---------------------------------------------------------------------------

def _call_chromosome_likelihood(args):
    """Worker: call ancestral states for one chromosome using likelihood."""
    (
        chrom, species_paths, tree_text,
        model_name, model_kappa, model_base_freqs, model_rates,
        high_thresh, low_thresh, out_path,
    ) = args

    tree = parse_newick(tree_text)
    model = build_model(
        model_name, kappa=model_kappa,
        base_freqs=model_base_freqs, rates=model_rates,
    )
    P_cache = _build_P_cache(tree, model)

    species_seqs: Dict[str, str] = {}
    length: Optional[int] = None
    for name, path in species_paths.items():
        _, seq = read_fasta(path)
        species_seqs[name] = seq
        if length is None:
            length = len(seq)
        elif len(seq) != length:
            raise ValueError(
                f"Length mismatch on {chrom}: expected {length}, got {len(seq)}"
            )

    anc = []
    for i in range(length):
        bases = {name: seq[i].upper() for name, seq in species_seqs.items()}
        anc.append(
            call_ancestral_base_likelihood(
                tree, bases, model, high_thresh, low_thresh, _P_cache=P_cache,
            )
        )

    write_fasta(out_path, f">{chrom}", "".join(anc))
    return chrom
