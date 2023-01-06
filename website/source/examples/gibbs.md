---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Gibbs free energy

```{code-cell} ipython3
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy import PhonopyQHA
from phonopy.units import CastepToTHz

import profess
import ase_tools
```

## Introduction

Lithium adopts the BCC structure in ambient conditions, but it favors the FCC structure at low temperatures. See the following paper for details:

> Ackland, et al., "Quantum and isotope effects in lithium metal," _Science_ **356**, 1254 (2017).

In this example, we compute the Gibbs free energy for the two structures across a range of temperatures, eventually predicting the phase transition.

+++

## Static lattice approximation

First, we define a helper function that creates a ``profess`` system with prescribed cell vectors and ion coordinates.

```{code-cell} ipython3
def build_system(box_vecs, xyz_coords):

    system = (
        profess.System.create(box_vecs, 1200, ['a','ev'])
        .add_ions('potentials/li.gga.recpot', xyz_coords, 'a')
        .add_electrons()
        .add_wang_teter_functional()              # kinetic energy
        .add_hartree_functional()                 # electrostatic energies
        .add_ion_electron_functional()
        .add_ion_ion_interaction()
        .add_perdew_burke_ernzerhof_functional()  # exchange-correlation energy
    )
    return system.minimize_energy()
```

Then, we find the equilibrium geometries for BCC Li and FCC Li, observing that BCC is roughly 1 meV/atom higher in energy than FCC. These calculations assume T=0 and a static lattice.

```{code-cell} ipython3
# initial guess for volume [A^3]
vol_per_atom = 20

# determine bcc energy
box_vec = (2*vol_per_atom)**(1.0/3.0) * np.array([[-0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5]])
system = build_system(box_vec, [(0,0,0)])
print('----------------------\nRelaxing BCC structure\n----------------------')
ase_tools.minimize_forces_stress(system, 'BFGSLineSearch', 1e-4)
print('\n    Final volume: {:4.2f}  A3'.format(system.volume('a3')))
print('    Final energy: {:4.3f} eV\n'.format(system.energy('ev')))
ene_bcc = system.energy('ev')

# determine fcc energy
box_vec = (4*vol_per_atom)**(1.0/3.0) * np.array([[0.5,0.5,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5]])
system = build_system(box_vec, [(0,0,0)])
print('----------------------\nRelaxing FCC structure\n----------------------')
ase_tools.minimize_forces_stress(system, 'BFGSLineSearch', 1e-4)
print('\n    Final volume: {:4.2f}  A3'.format(system.volume('a3')))
print('    Final energy: {:4.3f} eV\n'.format(system.energy('ev')))
ene_fcc = system.energy('ev')
```

## Gibbs free energy with Phonopy

Finally, we use ``Phonopy`` to compute Gibbs free energies with the quasi-harmonic approximation.
See the [Phonopy documentation](https://phonopy.github.io/phonopy/phonopy-module.html) for more details.

The thermal analysis explains why the bcc structure is observed at room temperature, despite fcc having slighly lower energy as $T\to 0$. The difference between the energies obtained in the previous section and the $T\to 0$ energies in the plot are the zero-point energies.

```{code-cell} ipython3
vols_per_atom = np.linspace(20,23,6)
temperatures = np.linspace(0,300,100)

for structure in ['bcc', 'fcc']:

    # create empty arrays
    volumes = np.zeros(vols_per_atom.size)
    energies = np.zeros(vols_per_atom.size)
    free_energies = np.zeros([temperatures.size, vols_per_atom.size])
    entropies = np.zeros([temperatures.size, vols_per_atom.size])
    heat_capacities = np.zeros([temperatures.size, vols_per_atom.size])

    for i, vol_per_atom in enumerate(vols_per_atom):

        # set box vectors
        if structure=='bcc':
            box_vec = (2*vol_per_atom)**(1.0/3.0) * np.array([[-0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5]])
        elif structure=='fcc':
            box_vec = (4*vol_per_atom)**(1.0/3.0) * np.array([[0.5,0.5,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5]])

        # compute volume and energy
        system = build_system(box_vec, [(0,0,0)])
        volumes[i] = system.volume('a3')
        energies[i] = system.energy('ev')

        # generate supercells with phonopy
        unitcell = PhonopyAtoms(symbols=(['Li']*1), cell=box_vec, scaled_positions=[(0,0,0)])
        supercell_matrix = [[2,0,0], [0,2,0], [0,0,2]]
        phonon = Phonopy(unitcell, supercell_matrix, primitive_matrix=np.identity(3), factor=CastepToTHz)
        phonon.generate_displacements(distance=0.01)
        disps = phonon.displacements
        supercells = phonon.supercells_with_displacements

        # compute forces for all supercells
        set_of_forces = []
        for supercell in supercells:
            system = build_system(supercell.cell, supercell.positions)
            forces = np.array(system.forces())*27.2114/0.529177
            drift_force = forces.sum(axis=0)
            for j, _ in enumerate(forces):
                forces[j] -= drift_force / forces.shape[0]
            set_of_forces.append(forces)

        # phonopy post-process
        phonon.produce_force_constants(forces=set_of_forces)

        # get thermal properties
        phonon.run_mesh([48,48,48])
        phonon.run_thermal_properties(temperatures=temperatures)
        free_energies[:,i] = phonon.get_thermal_properties_dict()['free_energy']
        entropies[:,i] = phonon.get_thermal_properties_dict()['entropy']
        heat_capacities[:,i] = phonon.get_thermal_properties_dict()['heat_capacity']
        
    # conduct quasi-harmonic analysis
    qha = PhonopyQHA(volumes=volumes,
                     electronic_energies=energies,
                     temperatures=temperatures,
                     free_energy=free_energies,
                     cv=heat_capacities,
                     entropy=entropies,
                     eos='vinet',
                     verbose=False)

    # extract gibbs energies
    if structure=='bcc':
        gibbs_bcc = qha.gibbs_temperature
    elif structure=='fcc':
        gibbs_fcc = qha.gibbs_temperature

plt.plot(temperatures[:-1], gibbs_bcc, label='bcc')
plt.plot(temperatures[:-1], gibbs_fcc, label='fcc')
plt.xlabel('Temperature [K]')
plt.ylabel('Gibbs Free Energy [eV/atom]')
plt.legend()
plt.show()
```
