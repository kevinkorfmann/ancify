ancify
======

**Ancestral allele polarization for any species.**

ancify is a config-driven Python pipeline that determines the ancestral state
at every position in a reference genome by comparing pairwise alignments from
multiple outgroup species. It uses a two-tier inner/outer outgroup voting
scheme with case-encoded confidence levels.

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

Get started in five minutes:

.. code-block:: bash

   pip install .
   ancify init -o config.yaml
   # edit config.yaml with your species, alignments, and paths
   ancify run -c config.yaml


Why ancify?
-----------

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
