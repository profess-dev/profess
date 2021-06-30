***************
Modify a system
***************

.. currentmodule:: profess.System

For fixed ion positions, one finds the electronic ground
state by varying the density to minimize the total energy.

.. autosummary::
    :toctree: generated
    :nosignatures:

    minimize_energy
    minimize_energy_tpsd

Other common modifications include moving the ions or changing the simulation box.

.. autosummary::
    :toctree: generated
    :nosignatures:

    move_ions
    set_box

Some uncommon modifications are useful for testing purposes.

.. autosummary::
    :toctree: generated
    :nosignatures:

    remove_functional

