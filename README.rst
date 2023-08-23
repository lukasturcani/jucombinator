Create all possible substitutions of a set of functional groups on a skeleton
structure.

Usage
=====

.. code-block:: python

  import jucombinator

  combinations = jucombinator.substitute(
      skeleton="c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1",
      substituents=["N(C)C", "O", "N", "S", "C", "F"],
      n=2,
  )

Intallation
===========

.. code-block:: bash

  pip install jucombinator
