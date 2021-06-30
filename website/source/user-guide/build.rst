**************
Build a system
**************

.. currentmodule:: profess.System

Create an empty system
======================

The first step is to create an empty system, ordinarily using :py:meth:`create`. ::

    system = profess.System.create(box_vectors, energy_cutoff, ['a','ev'])

The arguments provided to :py:meth:`create` are the lattice vectors (in the rows of an array-like object or list), the plane wave cutoff energy, and the units for the previous variables (given as a list).

.. autosummary::
    :toctree: generated
    :nosignatures:

    create
    create_from_grid_shape

.. _add-ions:

Add ions
========

The standard method is :py:meth:`add_ions`, which uses a pseudopotential
file. A method for adding bare Coulomb ions is also available.

.. autosummary::
    :toctree: generated
    :nosignatures:

    add_ions
    add_coulomb_ions

.. _add-electrons:

Add electrons
=============

.. autosummary::
    :toctree: generated
    :nosignatures:

    add_electrons

.. _add-energy-functionals:

Add energy functionals
======================

Kinetic
-------

.. autosummary::
    :toctree: generated
    :nosignatures:

    add_huang_carter_functional
    add_luo_karasiev_trickey_functional
    add_generic_nonlocal_a_functional
    add_perrot_functional
    add_smargiassi_madden_functional
    add_thomas_fermi_functional
    add_wang_teter_functional
    add_wang_govind_carter_functional
    add_wang_govind_carter_1999_i_functional
    add_weizsaecker_functional

Electrostatic
-------------

.. autosummary::
    :toctree: generated
    :nosignatures:

    add_hartree_functional
    add_ion_electron_functional
    add_ion_ion_interaction

Exchange-correlation
--------------------

.. autosummary::
    :toctree: generated
    :nosignatures:

    add_libxc_functional
    add_perdew_burke_ernzerhof_functional
    add_perdew_zunger_functional

System complete
===============

After building the system, the usual next step is to find the ground state
with :py:meth:`minimize_energy`.

One may then :doc:`inspect the system <inspect>`, computing or reporting important properties, or
:doc:`modify the system <modify>`, such as by moving the ions.
