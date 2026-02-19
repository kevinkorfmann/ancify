ancify
======

**Infer ancestral alleles for any species using outgroup alignments.**

ancify is a config-driven Python pipeline that determines the ancestral state
at every position in a reference genome by comparing pairwise alignments from
multiple outgroup species.  It uses a two-tier inner/outer outgroup voting
scheme with case-encoded confidence levels.

.. code-block:: bash

   pip install .
   ancify init -o config.yaml
   # edit config.yaml
   ancify run -c config.yaml

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   quickstart
   installation
   configuration
   cli

.. toctree::
   :maxdepth: 2
   :caption: Concepts

   algorithm
   species_guide
   evaluation

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   changelog


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
