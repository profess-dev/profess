import numpy as np
import os
import sys

sys.path.append(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), '../../../../external/ase/'))
from ase import Atoms
import ase.io

import profess
import ase_tools

to_eV_per_A3 = 1.0 / 160.2177 # from GPa
to_Ha_per_B3 = 1.0 / 29421.02648 # from GPa

#------------------------------------------------

# Collect command line input and load random structure
pressure = float(sys.argv[1])
seed = sys.argv[2]
element = seed.split('-')[0][2:]
atoms = ase.io.read(seed + '.xyz')
box_vectors = atoms.cell[...]
xyz_coords = atoms.get_positions()

# Relax in stages to refresh the grid periodically
stage = 0
num_converged = 0
while stage < 100:
    stage += 1
    print('---------- STAGE {} ----------'.format(stage))
    steps_per_stage = 20
    if stage <= 3:  # use short stages at first
        steps_per_stage = 3

    # Build system and relax geometry
    energy_cutoff = 600
    potential = '../potentials/' + element.lower() + '.gga.recpot'
    system = (
        profess.System.create(box_vectors, energy_cutoff, ['a','ev'])
        .add_ions(potential, xyz_coords, 'a')
        .add_electrons()
        .add_generic_nonlocal_a_functional(       # kinetic energy
            (5+np.sqrt(5))/6,
            (5-np.sqrt(5))/6,
            lambda x: np.exp(x),
            lambda x: np.exp(x))
        .add_hartree_functional()                 # electrostatic energies
        .add_ion_electron_functional()
        .add_ion_ion_interaction()
        .add_perdew_burke_ernzerhof_functional()  # exchange-correlation
    ).minimize_energy()
    converged = ase_tools.minimize_forces_stress(
        system,
        'BFGSLineSearch',
        tol=2e-4,
        steps=steps_per_stage,
        scalar_pressure=pressure*to_Ha_per_B3)
    box_vectors = np.array(system.box_vectors('a'))
    xyz_coords = np.array(system.ions_xyz_coords('a'))

    # If converged twice, output final structure and fake castep file
    if converged:
        num_converged += 1
    if num_converged == 2:
        atoms = Atoms(element+str(xyz_coords.shape[0]),
                      positions=xyz_coords,
                      cell=box_vectors),
        ase.io.write(seed + '-out.xyz', atoms, 'extxyz')
        volume = system.volume('a3')
        pressure = system.pressure('gpa')
        enthalpy = system.enthalpy('ev')
        with open(seed + '.castep', 'w') as f:
            f.write("Current cell volume = {:25.15f} A**3\n".format(volume))
            f.write("*  Pressure:   {:25.15f}\n".format(pressure))
            f.write("Python: Final Enthalpy     = {:25.15f} eV\n".format(enthalpy))
        break
