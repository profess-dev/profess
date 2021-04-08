**************
Build a system
**************

.. currentmodule:: profess.System

Systems are built in the following sequence:

* :ref:`create the system <create-the-system>`;
* :ref:`add ions <add-ions>`;
* :ref:`add electrons <add-electrons>`;
* :ref:`add energy functionals <add-energy-functionals>`.

.. _create-the-system:

Create the system
=================

Systems are created as follows. ::

    shape = profess.System.get_shape(box_vectors, energy_cutoff, ('a','ev'))
    system = profess.System(shape)
    system.set_box(box_vectors, 'a')

.. _add-ions:

Add ions
========

The standard method for adding ions to a system is :py:meth:`add_ions`, which uses a pseudopotential
file. But a method for adding bare Coulomb ions is also available.

.. autosummary::
    :toctree: generated
    :nosignatures:

    add_ions
    add_coulomb_ions
    add_harmonic_ions

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

.. autosummary::
    :toctree: generated
    :nosignatures:

    add_hartree_functional
    add_ion_electron_functional
    add_luo_karasiev_trickey_functional
    add_libxc_functional
    add_kinetic_class_a_functional
    add_perdew_burke_ernzerhof_functional
    add_perdew_zunger_functional
    add_perrot_functional
    add_smargiassi_madden_functional
    add_thomas_fermi_functional
    add_wang_teter_functional
    add_wang_govind_carter_1999_i_functional
    add_weizsaecker_functional
