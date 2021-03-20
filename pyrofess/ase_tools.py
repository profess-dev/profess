# Copyright (c) 2019-2020 William C. Witt
# SPDX-License-Identifier: MIT

import numpy as np
import os
import sys

from ase import Atoms
from ase.optimize import BFGS, BFGSLineSearch, LBFGS, LBFGSLineSearch
from ase.calculators.profess import Profess
from ase.constraints import FixAtoms, ExpCellFilter

def init_optimizer(atoms, algorithm):

    if algorithm == 'BFGS':
        optimizer = BFGS(atoms)
    elif algorithm == 'BFGSLineSearch':
        optimizer = BFGSLineSearch(atoms)
    elif algorithm == 'LBFGS':
        optimizer = LBFGS(atoms)
    elif algorithm == 'LBFGSLineSearch':
        optimizer = LBFGSLineSearch(atoms)
    else:
        sys.exit('algorithm not recognized')
    return optimizer

def minimize_forces(system, algorithm, fmax=1e-4, steps=100):
    atoms = Atoms('H' + str(system.ions.count()),
            positions=system.ions.xyz_coords(),
            cell=system.box_vectors(),
            pbc=True)
    atoms.calc = Profess(system)
    optimizer = init_optimizer(atoms, algorithm)
    optimizer.run(fmax, steps)
    return optimizer.converged()
    
def minimize_stress(
        system,
        algorithm,
        tol=1e-4,
        steps=100,
        hydrostatic_strain=False,
        constant_volume=False,
        scalar_pressure=0.0):
    atoms = Atoms('H' + str(system.ions.count()),
            positions=system.ions.xyz_coords(),
            cell=system.box_vectors(),
            pbc=True)
    atoms.calc = Profess(system)
    atoms.set_constraint(FixAtoms(mask=[True for atom in atoms]))
    ecf = ExpCellFilter(
            atoms,
            hydrostatic_strain=hydrostatic_strain,
            constant_volume=constant_volume,
            scalar_pressure=scalar_pressure)
    optimizer = init_optimizer(ecf, algorithm)
    optimizer.run(tol, steps)
    return optimizer.converged()

def minimize_forces_stress(
        system,
        algorithm,
        tol=1e-4,
        steps=100,
        hydrostatic_strain=False,
        constant_volume=False,
        scalar_pressure=0.0):
    atoms = Atoms('H' + str(system.ions.count()),
            positions=system.ions_xyz_coords(),
            cell=system.box_vectors(),
            pbc=True)
    atoms.calc = Profess(system)
    ecf = ExpCellFilter(
            atoms,
            hydrostatic_strain=hydrostatic_strain,
            constant_volume=constant_volume,
            scalar_pressure=scalar_pressure)
    optimizer = init_optimizer(ecf, algorithm)
    optimizer.run(tol, steps)
    return optimizer.converged()

def optimize_geometry(init_system, box_vectors, xyz_coords, energy_cutoff):
    converged = False
    stages = 0
    while not converged:
        stages += 1
        if stages <= 10:
            steps_per_stage = 2
        else:
            steps_per_stage = 20
        print('steps_per_stage is', steps_per_stage)
        system = init_system(box_vectors.T, xyz_coords, energy_cutoff)
        converged = minimize_forces_stress(
                system, 'BFGSLineSearch', steps=steps_per_stage)
        print('')
        print('stages completed:', stages)
        print('')
        if converged:
            print('converged! final energy:', system.energy())
            return system
        else:
            box_vectors = np.array(system.box.vectors()).T
            xyz_coords = np.array(system.ions.xyz_coords())
