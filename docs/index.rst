
.. container:: epigraph

   Was einst im tiefen Schweigen schlief,
   ist nun ein fernes Echo seiner selbst.

**Ancestral allele polarization for any species.**

ancify is a config-driven Python pipeline that determines the ancestral state
at every position in a reference genome by comparing pairwise alignments from
multiple outgroup species. It supports four inference methods: **two-tier voting**,
**Fitch parsimony** on a phylogenetic tree, **likelihood** (Felsenstein pruning with
substitution models), and a **machine-learning classifier** (LightGBM) — all
producing case-encoded confidence levels.

.. code-block:: text

                      ┌──────────┐
       Inner:         │  Bonobo  │──┐
       closely        ├──────────┤  │  majority   ┌──────────────┐
       related        │  Chimp   │──┼────vote────▶│    Inner     │
       species        ├──────────┤  │             │  consensus   │
                      │ Gorilla  │──┘             └──────┬───────┘
                                                         │ compare
                      ┌──────────┐              ┌────────▼───────┐
       Outer:         │ Macaque  │──────vote────▶│   Ancestral   │
       distant        └──────────┘              │     call      │
       species                                  └───────────────┘

                         Agree? → UPPERCASE (high confidence)
                      One tier? → lowercase (low confidence)
                      Disagree? → n (unresolved)
                    Both empty? → N (missing)

Four methods, one config field:

.. code-block:: yaml

   method: voting      # default — two-tier majority vote (above diagram)
   method: parsimony   # Fitch algorithm on a Newick phylogenetic tree
   method: likelihood  # Felsenstein pruning + JC69/K80/HKY85/GTR (tree with branch lengths)
   method: ml          # LightGBM classifier trained on your alignment data

Get started in five minutes:

.. code-block:: bash

   pip install .
   ancify init -o config.yaml
   # edit config.yaml with your species, alignments, and paths
   ancify run -c config.yaml


Why ancify?
-----------

- **Four inference methods** — choose the approach that fits your data: two-tier voting (fast, transparent), Fitch parsimony (tree-aware, resolves ILS), likelihood (substitution models + branch lengths, posterior probabilities), or ML classifier (learns substitution biases from your data). See :doc:`algorithm`.
- **GPU-accelerated** — optional PyTorch backend turns 12-hour genome-wide runs into ~2 minutes on NVIDIA GPUs. See :doc:`performance`.
- **Species-agnostic** — works with humans, mice, flies, fish, plants, or any species with outgroup alignments.
- **Educational** — the docs teach you the population genetics *behind* polarization, not just the buttons to press.
- **Config-driven** — one YAML file controls everything. No scripts to edit.
- **Transparent** — confidence is encoded directly in the output (uppercase/lowercase). You always know how reliable each call is.
- **Validated** — tested against the Ensembl EPO 13-primate ancestral reference with >99.9% agreement.

.. toctree::
   :maxdepth: 2
   :caption: Learn

   background
   quickstart
   tutorials

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   configuration
   performance
   cli
   species_guide

.. toctree::
   :maxdepth: 2
   :caption: Deep Dives

   algorithm
   evaluation
   faq

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   glossary
   changelog


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
