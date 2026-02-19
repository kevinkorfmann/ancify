"""Fitch parsimony for ancestral allele reconstruction on a phylogenetic tree.

Implements the Fitch (1971) algorithm:

1. **Bottom-up (post-order):** At each leaf, assign the observed allele as a
   singleton set.  ``N`` is treated as ``{A, C, G, T}`` (compatible with
   everything).  At each internal node, take the intersection of children's
   sets if non-empty, otherwise their union.

2. **Top-down (pre-order):** Starting at the root, pick a concrete allele from
   the node's set (preferring the parent's assignment for determinism, breaking
   ties alphabetically).  Propagate down.

The root allele is the inferred ancestral state.  Ambiguity at the root
(set size > 1 after the bottom-up pass) indicates multiple equally
parsimonious reconstructions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

VALID_ALLELES_SET = frozenset({"A", "C", "G", "T"})
VALID_BASES = VALID_ALLELES_SET


# ---------------------------------------------------------------------------
# Tree data structure
# ---------------------------------------------------------------------------

@dataclass
class TreeNode:
    """A node in an unrooted/rooted phylogenetic tree.

    Leaves have a *name* matching an outgroup species identifier.
    Internal nodes have ``name=None``.  Branch lengths are parsed but
    not used by the Fitch algorithm.
    """

    name: Optional[str] = None
    children: List["TreeNode"] = field(default_factory=list)
    branch_length: Optional[float] = None

    @property
    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def leaf_names(self) -> List[str]:
        """Return all leaf names in pre-order."""
        if self.is_leaf:
            return [self.name] if self.name else []
        names: List[str] = []
        for child in self.children:
            names.extend(child.leaf_names())
        return names


def get_leaf_names(tree: TreeNode) -> List[str]:
    """Return all leaf names from *tree* (module-level convenience wrapper)."""
    return tree.leaf_names()


# ---------------------------------------------------------------------------
# Newick parser
# ---------------------------------------------------------------------------

def parse_newick(text: str) -> TreeNode:
    """Parse a Newick-format string into a :class:`TreeNode` tree.

    Supports optional branch lengths (``name:length``) and nested clades.
    Whitespace is ignored.  The trailing semicolon is optional.

    Examples::

        >>> tree = parse_newick("((A,B),C);")
        >>> tree.leaf_names()
        ['A', 'B', 'C']
    """
    text = text.strip().rstrip(";").strip()
    if not text:
        raise ValueError("Empty Newick string")
    node, pos = _parse_node(text, 0)
    remaining = text[pos:].strip()
    if remaining:
        raise ValueError(
            f"Unexpected trailing characters in Newick string: {remaining!r}"
        )
    return node


def _parse_node(text: str, pos: int) -> Tuple[TreeNode, int]:
    """Recursive descent parser for a single Newick node starting at *pos*."""
    node = TreeNode()

    if pos < len(text) and text[pos] == "(":
        pos += 1  # consume '('
        while True:
            pos = _skip_whitespace(text, pos)
            child, pos = _parse_node(text, pos)
            node.children.append(child)
            pos = _skip_whitespace(text, pos)
            if pos < len(text) and text[pos] == ",":
                pos += 1  # consume ','
            else:
                break
        if pos >= len(text) or text[pos] != ")":
            raise ValueError(
                f"Expected ')' at position {pos} in Newick string"
            )
        pos += 1  # consume ')'

    pos = _skip_whitespace(text, pos)
    name, pos = _parse_label(text, pos)
    if name:
        node.name = name

    pos = _skip_whitespace(text, pos)
    if pos < len(text) and text[pos] == ":":
        pos += 1  # consume ':'
        length_str, pos = _parse_number(text, pos)
        node.branch_length = float(length_str) if length_str else None

    return node, pos


def _skip_whitespace(text: str, pos: int) -> int:
    while pos < len(text) and text[pos] in (" ", "\t", "\n", "\r"):
        pos += 1
    return pos


def _parse_label(text: str, pos: int) -> Tuple[str, int]:
    """Parse an unquoted or single-quoted Newick label."""
    if pos < len(text) and text[pos] == "'":
        pos += 1
        start = pos
        while pos < len(text) and text[pos] != "'":
            pos += 1
        label = text[start:pos]
        if pos < len(text):
            pos += 1  # consume closing quote
        return label, pos

    start = pos
    stop_chars = {"(", ")", ",", ":", ";", " ", "\t", "\n", "\r"}
    while pos < len(text) and text[pos] not in stop_chars:
        pos += 1
    return text[start:pos], pos


def _parse_number(text: str, pos: int) -> Tuple[str, int]:
    """Parse a numeric literal (branch length)."""
    pos = _skip_whitespace(text, pos)
    start = pos
    while pos < len(text) and text[pos] in "0123456789.eE+-":
        pos += 1
    return text[start:pos], pos


# ---------------------------------------------------------------------------
# Fitch algorithm
# ---------------------------------------------------------------------------

def fitch_bottom_up(
    node: TreeNode,
    leaf_alleles: Dict[str, str],
) -> Dict[int, Set[str]]:
    """Post-order traversal: compute the Fitch set at every node.

    Parameters
    ----------
    node : TreeNode
        Root of the (sub)tree.
    leaf_alleles : dict
        Mapping of leaf name to observed allele (single uppercase character).
        ``'N'`` or any character not in ``{A, C, G, T}`` is treated as
        the full set (wildcard).

    Returns
    -------
    dict
        Mapping of ``id(node)`` to a frozenset of possible alleles.
    """
    sets: Dict[int, Set[str]] = {}
    _fitch_up(node, leaf_alleles, sets)
    return sets


def _fitch_up(
    node: TreeNode,
    leaf_alleles: Dict[str, str],
    sets: Dict[int, Set[str]],
) -> Set[str]:
    if node.is_leaf:
        allele = leaf_alleles.get(node.name, "N").upper()
        if allele in VALID_ALLELES_SET:
            s = {allele}
        else:
            s = set(VALID_ALLELES_SET)
        sets[id(node)] = s
        return s

    child_sets = [_fitch_up(c, leaf_alleles, sets) for c in node.children]

    intersection = child_sets[0]
    for cs in child_sets[1:]:
        intersection = intersection & cs

    if intersection:
        sets[id(node)] = intersection
        return intersection
    else:
        union = set()
        for cs in child_sets:
            union = union | cs
        sets[id(node)] = union
        return union


def fitch_top_down(
    node: TreeNode,
    node_sets: Dict[int, Set[str]],
    parent_allele: Optional[str] = None,
) -> Dict[int, str]:
    """Pre-order traversal: assign a concrete allele at each node.

    Deterministic tie-breaking: prefer the parent's allele if it is in
    the node's Fitch set; otherwise pick the lexicographically smallest.

    Returns
    -------
    dict
        Mapping of ``id(node)`` to the assigned allele character.
    """
    assignments: Dict[int, str] = {}
    _fitch_down(node, node_sets, parent_allele, assignments)
    return assignments


def _fitch_down(
    node: TreeNode,
    node_sets: Dict[int, Set[str]],
    parent_allele: Optional[str],
    assignments: Dict[int, str],
) -> None:
    s = node_sets[id(node)]
    if parent_allele is not None and parent_allele in s:
        chosen = parent_allele
    else:
        chosen = min(sorted(s))
    assignments[id(node)] = chosen

    for child in node.children:
        _fitch_down(child, node_sets, chosen, assignments)


def fitch_ancestral(
    tree: TreeNode,
    leaf_alleles: Dict[str, str],
) -> Tuple[str, bool]:
    """Run the full Fitch algorithm and return the inferred root state.

    Parameters
    ----------
    tree : TreeNode
        The phylogenetic tree of outgroup species.
    leaf_alleles : dict
        Mapping of leaf name to observed allele at this position.

    Returns
    -------
    tuple of (str, bool)
        ``(root_allele, is_ambiguous)`` where *is_ambiguous* is ``True``
        when the root's Fitch set contained more than one allele.
    """
    node_sets = fitch_bottom_up(tree, leaf_alleles)
    root_set = node_sets[id(tree)]

    if root_set == VALID_ALLELES_SET:
        return "N", True

    assignments = fitch_top_down(tree, node_sets)
    root_allele = assignments[id(tree)]
    is_ambiguous = len(root_set) > 1

    return root_allele, is_ambiguous
