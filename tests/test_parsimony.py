"""Tests for ancify.parsimony -- Newick parsing and Fitch algorithm."""

import pytest

from ancify.parsimony import (
    TreeNode,
    fitch_ancestral,
    fitch_bottom_up,
    fitch_top_down,
    get_leaf_names,
    parse_newick,
)


# ======================================================================
# Newick parsing
# ======================================================================

class TestParseNewick:
    """Test the recursive-descent Newick parser."""

    def test_two_leaves(self):
        tree = parse_newick("(A,B);")
        assert sorted(tree.leaf_names()) == ["A", "B"]
        assert len(tree.children) == 2

    def test_three_leaves(self):
        tree = parse_newick("((A,B),C);")
        names = tree.leaf_names()
        assert sorted(names) == ["A", "B", "C"]
        assert len(tree.children) == 2
        assert len(tree.children[0].children) == 2  # (A,B) clade

    def test_four_species_primate_tree(self):
        tree = parse_newick("(((bonobo,chimp),gorilla),macaque);")
        assert sorted(tree.leaf_names()) == [
            "bonobo", "chimp", "gorilla", "macaque"
        ]

    def test_branch_lengths_parsed(self):
        tree = parse_newick("(A:0.1,B:0.2):0.05;")
        assert tree.branch_length == pytest.approx(0.05)
        assert tree.children[0].branch_length == pytest.approx(0.1)
        assert tree.children[1].branch_length == pytest.approx(0.2)

    def test_branch_lengths_ignored_by_fitch(self):
        tree_with = parse_newick("(A:0.1,B:0.2);")
        tree_without = parse_newick("(A,B);")
        alleles = {"A": "G", "B": "G"}
        assert fitch_ancestral(tree_with, alleles) == fitch_ancestral(tree_without, alleles)

    def test_no_semicolon(self):
        tree = parse_newick("(A,B)")
        assert sorted(tree.leaf_names()) == ["A", "B"]

    def test_whitespace_tolerance(self):
        tree = parse_newick("  ( A , B )  ;  ")
        assert sorted(tree.leaf_names()) == ["A", "B"]

    def test_single_leaf(self):
        tree = parse_newick("A;")
        assert tree.leaf_names() == ["A"]
        assert tree.is_leaf

    def test_quoted_labels(self):
        tree = parse_newick("('species one','species two');")
        assert sorted(tree.leaf_names()) == ["species one", "species two"]

    def test_empty_string_raises(self):
        with pytest.raises(ValueError, match="Empty"):
            parse_newick("")

    def test_unbalanced_parentheses_raises(self):
        with pytest.raises(ValueError):
            parse_newick("((A,B)")

    def test_deeply_nested(self):
        tree = parse_newick("((((A,B),C),D),E);")
        assert sorted(tree.leaf_names()) == ["A", "B", "C", "D", "E"]

    def test_multifurcation(self):
        tree = parse_newick("(A,B,C,D);")
        assert sorted(tree.leaf_names()) == ["A", "B", "C", "D"]
        assert len(tree.children) == 4

    def test_internal_node_labels(self):
        tree = parse_newick("((A,B)ancestor1,C);")
        inner = tree.children[0]
        assert inner.name == "ancestor1"
        assert sorted(inner.leaf_names()) == ["A", "B"]


class TestGetLeafNames:
    """Test the module-level get_leaf_names convenience function."""

    def test_matches_method(self):
        tree = parse_newick("((A,B),C);")
        assert get_leaf_names(tree) == tree.leaf_names()

    def test_preserves_order(self):
        tree = parse_newick("((X,Y),Z);")
        assert get_leaf_names(tree) == ["X", "Y", "Z"]


# ======================================================================
# Fitch bottom-up pass
# ======================================================================

class TestFitchBottomUp:
    """Test the post-order (bottom-up) pass of the Fitch algorithm."""

    def test_all_same_allele(self):
        tree = parse_newick("(A,B);")
        sets = fitch_bottom_up(tree, {"A": "G", "B": "G"})
        root_set = sets[id(tree)]
        assert root_set == {"G"}

    def test_different_alleles_union(self):
        tree = parse_newick("(A,B);")
        sets = fitch_bottom_up(tree, {"A": "G", "B": "T"})
        root_set = sets[id(tree)]
        assert root_set == {"G", "T"}

    def test_intersection_propagates(self):
        tree = parse_newick("((A,B),C);")
        sets = fitch_bottom_up(tree, {"A": "G", "B": "G", "C": "G"})
        root_set = sets[id(tree)]
        assert root_set == {"G"}

    def test_missing_data_is_wildcard(self):
        tree = parse_newick("(A,B);")
        sets = fitch_bottom_up(tree, {"A": "G", "B": "N"})
        root_set = sets[id(tree)]
        assert root_set == {"G"}

    def test_all_missing_is_full_set(self):
        tree = parse_newick("(A,B);")
        sets = fitch_bottom_up(tree, {"A": "N", "B": "N"})
        root_set = sets[id(tree)]
        assert root_set == {"A", "C", "G", "T"}

    def test_primate_tree_all_agree(self):
        tree = parse_newick("(((bonobo,chimp),gorilla),macaque);")
        alleles = {
            "bonobo": "A", "chimp": "A", "gorilla": "A", "macaque": "A"
        }
        sets = fitch_bottom_up(tree, alleles)
        assert sets[id(tree)] == {"A"}

    def test_primate_tree_ils_scenario(self):
        """Bonobo+chimp share G, gorilla+macaque share A.

        Tree: (((bonobo,chimp),gorilla),macaque)
        Bottom-up:
          (bonobo,chimp) -> {G} (intersection)
          ((bonobo,chimp),gorilla) -> {G} | {A} = {A,G} (union)
          root -> {A,G} & {A} = {A} (intersection with macaque)
        """
        tree = parse_newick("(((bonobo,chimp),gorilla),macaque);")
        alleles = {
            "bonobo": "G", "chimp": "G", "gorilla": "A", "macaque": "A"
        }
        sets = fitch_bottom_up(tree, alleles)
        assert sets[id(tree)] == {"A"}


# ======================================================================
# Fitch top-down pass
# ======================================================================

class TestFitchTopDown:
    """Test the pre-order (top-down) pass of the Fitch algorithm."""

    def test_unambiguous_propagation(self):
        tree = parse_newick("(A,B);")
        sets = fitch_bottom_up(tree, {"A": "C", "B": "C"})
        assignments = fitch_top_down(tree, sets)
        assert assignments[id(tree)] == "C"

    def test_deterministic_tiebreak(self):
        tree = parse_newick("(A,B);")
        sets = fitch_bottom_up(tree, {"A": "G", "B": "T"})
        assignments = fitch_top_down(tree, sets)
        assert assignments[id(tree)] == "G"

    def test_parent_preference(self):
        """If parent allele is in child's set, child inherits it."""
        tree = parse_newick("((A,B),C);")
        alleles = {"A": "G", "B": "T", "C": "G"}
        sets = fitch_bottom_up(tree, alleles)
        assignments = fitch_top_down(tree, sets)
        root_allele = assignments[id(tree)]
        assert root_allele == "G"


# ======================================================================
# Full Fitch algorithm (fitch_ancestral)
# ======================================================================

class TestFitchAncestral:
    """Test the complete two-pass Fitch algorithm."""

    def test_all_agree_unambiguous(self):
        tree = parse_newick("(A,B,C);")
        allele, ambiguous = fitch_ancestral(tree, {"A": "T", "B": "T", "C": "T"})
        assert allele == "T"
        assert ambiguous is False

    def test_two_vs_one_unambiguous(self):
        tree = parse_newick("((A,B),C);")
        allele, ambiguous = fitch_ancestral(
            tree, {"A": "A", "B": "A", "C": "G"}
        )
        assert ambiguous is True

    def test_star_topology_two_vs_one(self):
        """On a star tree (A,B,C), 2×A + 1×G is ambiguous at the root.

        Bottom-up: {A} ∩ {A} ∩ {G} = {} → union {A,G}.
        Without internal structure the tree cannot resolve 2-vs-1.
        """
        tree = parse_newick("(A,B,C);")
        allele, ambiguous = fitch_ancestral(
            tree, {"A": "A", "B": "A", "C": "G"}
        )
        assert allele in ("A", "G")
        assert ambiguous is True

    def test_primate_ils_resolves_correctly(self):
        """The classic example: tree structure resolves what flat voting cannot.

        Tree: (((bonobo,chimp),gorilla),macaque)
        bonobo=G, chimp=G, gorilla=A, macaque=A

        Parsimony: A is ancestral (1 mutation: G arose in bonobo-chimp ancestor)
        Voting would be ambiguous (2 vs 2).
        """
        tree = parse_newick("(((bonobo,chimp),gorilla),macaque);")
        alleles = {
            "bonobo": "G", "chimp": "G", "gorilla": "A", "macaque": "A"
        }
        allele, ambiguous = fitch_ancestral(tree, alleles)
        assert allele == "A"
        assert ambiguous is False

    def test_all_different_four_species(self):
        """All four leaves have different alleles: completely uninformative.

        Bottom-up: (A,B) → {A}∪{C}={A,C}; (C,D) → {G}∪{T}={G,T}
        root → {A,C} ∩ {G,T} = {} → {A,C,G,T} = full set → treated as N.
        """
        tree = parse_newick("((A,B),(C,D));")
        allele, ambiguous = fitch_ancestral(
            tree, {"A": "A", "B": "C", "C": "G", "D": "T"}
        )
        assert ambiguous is True
        assert allele == "N"

    def test_missing_one_leaf(self):
        tree = parse_newick("((A,B),C);")
        allele, ambiguous = fitch_ancestral(
            tree, {"A": "G", "B": "N", "C": "G"}
        )
        assert allele == "G"
        assert ambiguous is False

    def test_all_missing_returns_N(self):
        tree = parse_newick("(A,B,C);")
        allele, ambiguous = fitch_ancestral(
            tree, {"A": "N", "B": "N", "C": "N"}
        )
        assert allele == "N"
        assert ambiguous is True

    def test_single_species_tree(self):
        tree = parse_newick("A;")
        allele, ambiguous = fitch_ancestral(tree, {"A": "C"})
        assert allele == "C"
        assert ambiguous is False

    def test_single_species_missing(self):
        tree = parse_newick("A;")
        allele, ambiguous = fitch_ancestral(tree, {"A": "N"})
        assert allele == "N"
        assert ambiguous is True

    def test_lowercase_input_handled(self):
        tree = parse_newick("(A,B);")
        allele, ambiguous = fitch_ancestral(tree, {"A": "g", "B": "g"})
        assert allele == "G"
        assert ambiguous is False

    def test_missing_species_in_alleles_treated_as_N(self):
        tree = parse_newick("(A,B);")
        allele, ambiguous = fitch_ancestral(tree, {"A": "T"})
        assert allele == "T"
        assert ambiguous is False

    def test_five_species_deep_tree(self):
        tree = parse_newick("((((A,B),C),D),E);")
        alleles = {"A": "G", "B": "G", "C": "G", "D": "G", "E": "G"}
        allele, ambiguous = fitch_ancestral(tree, alleles)
        assert allele == "G"
        assert ambiguous is False

    def test_five_species_one_derived(self):
        tree = parse_newick("((((A,B),C),D),E);")
        alleles = {"A": "T", "B": "G", "C": "G", "D": "G", "E": "G"}
        allele, ambiguous = fitch_ancestral(tree, alleles)
        assert allele == "G"

    def test_symmetric_tree_two_vs_two(self):
        """On ((A,B),(C,D)): if A,B=G and C,D=T, root is ambiguous."""
        tree = parse_newick("((A,B),(C,D));")
        alleles = {"A": "G", "B": "G", "C": "T", "D": "T"}
        allele, ambiguous = fitch_ancestral(tree, alleles)
        assert ambiguous is True
        assert allele in ("G", "T")

    def test_asymmetric_tree_two_vs_two_resolves(self):
        """On (((A,B),C),D): if A,B=G, C=T, D=T, root set is {T}.

        Bottom-up:
          (A,B) -> {G}
          ((A,B),C) -> {G} | {T} = {G,T} (union)
          root -> {G,T} & {T} = {T} (intersection with D)
        """
        tree = parse_newick("(((A,B),C),D);")
        alleles = {"A": "G", "B": "G", "C": "T", "D": "T"}
        allele, ambiguous = fitch_ancestral(tree, alleles)
        assert allele == "T"
        assert ambiguous is False


# ======================================================================
# TreeNode properties
# ======================================================================

class TestTreeNode:
    """Test TreeNode dataclass behavior."""

    def test_leaf_detection(self):
        leaf = TreeNode(name="A")
        assert leaf.is_leaf is True

    def test_internal_detection(self):
        node = TreeNode(children=[TreeNode(name="A"), TreeNode(name="B")])
        assert node.is_leaf is False

    def test_leaf_names_empty_internal(self):
        node = TreeNode()
        assert node.leaf_names() == []

    def test_leaf_names_preserves_preorder(self):
        tree = parse_newick("((A,B),(C,D));")
        assert tree.leaf_names() == ["A", "B", "C", "D"]
