API Reference
=============

ancify can be used as a Python library for programmatic access to the
ancestral calling pipeline.

Quick example
-------------

.. code-block:: python

   from ancify.config import load_config
   from ancify.project import run_projection
   from ancify.ancestral import run_ancestral_calling, call_ancestral_base

   # Run the full pipeline
   cfg = load_config("config.yaml")
   run_projection(cfg)
   run_ancestral_calling(cfg)

   # Or call the core function directly
   base = call_ancestral_base(
       inner_bases=["A", "A", "G"],
       outer_bases=["A"],
   )
   # Returns "A" (high confidence)


ancify.utils
------------

.. automodule:: ancify.utils
   :members:
   :undoc-members:


ancify.config
-------------

.. automodule:: ancify.config
   :members:
   :undoc-members:


ancify.project
--------------

.. automodule:: ancify.project
   :members:
   :undoc-members:


ancify.ancestral
----------------

.. automodule:: ancify.ancestral
   :members:
   :undoc-members:


ancify.evaluate
---------------

.. automodule:: ancify.evaluate
   :members:
   :undoc-members:


ancify.cli
----------

.. automodule:: ancify.cli
   :members:
   :undoc-members:
