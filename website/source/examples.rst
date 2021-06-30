********
Examples
********

.. toctree::
    :caption: Basic
    :hidden:

    examples/aluminum
    examples/h-and-h2
    examples/magnesium

.. toctree::
    :caption: More advanced
    :hidden:

    examples/gibbs
    examples/airss/airss


Basic
=====

:doc:`Aluminum <examples/aluminum>`
    Extended version of the basic example shown on the
    :doc:`User Guide <user-guide>` page.

:doc:`H and H2 <examples/h-and-h2>`
    For these simple examples, the exact
    kinetic energy functional is known.

:doc:`Magnesium <examples/magnesium>`
    Obtains relaxed lattice parameters for hexagonal close-packed
    (HCP) magnesium using several methods.

More Advanced
=============

:doc:`Gibbs free energy <examples/gibbs>`
    Using Phonopy and the quasi-harmonic approximations, computes the
    Gibbs free energy over a range of temperatures for fcc and bcc lithium.
    Demonstrates why the bcc structure is observed at room temperature, despite
    having slighly higher energy when the lattice is assumed static.

:doc:`Ab initio random structure searching (AIRSS) <examples/airss/airss>`
    Generates random unit cells with four magnesium atoms and relaxes them,
    following the AIRSS protocol for crystal structure prediction.
