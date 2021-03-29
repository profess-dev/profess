# Copyright (c) 2019-2020 William C. Witt
# SPDX-License-Identifier: MIT

import numpy as np
import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), '../pyrofess/'))
import ase_tools
from tools_for_tests import murnaghan
sys.path.append(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), '../build/'))
import profess

class TestGeomOpt(unittest.TestCase):

    def test_bcc_lithium(self):

        # create system and compute ground state energy
        system = profess.System([15,15,15])
        box_vectors = 3.48*np.eye(3)
        system.set_box(box_vectors, 'a')
        system.add_ions(
            'potentials/li.gga.recpot',
            box_vectors[0,0]*np.array([[0,0,0],[0.5,0.5,0.5]]),
            'a')
        (
        system
               .add_hartree_functional()
               .add_ion_electron_functional()
               .add_perdew_burke_ernzerhof_functional()
               .add_wang_teter_functional()
        )
        system.add_electrons(system.total_ion_charge())
        system.minimize_energy()
        energy = system.energy()
        # peturb ions, then restore by minimizing forces, keeping box fixed
        system.move_ions(
            box_vectors[0,0]*np.array([[0.0,0.1,0.0],[0.6,0.4,0.6]]),
            'a')
        ase_tools.minimize_forces(system, 'BFGSLineSearch', 1e-4)
        self.assertAlmostEqual(energy, system.energy(), places=5)
        # predict relaxed energy by fitting to murnaghan equation
        lengths = 3.48*np.linspace(0.98,1.02,9)
        energies = []
        volumes = []
        for length in lengths:
            system.set_box(length*np.eye(3), 'a')
            system.minimize_energy()
            energies.append(system.energy() / 2.0)
            volumes.append(system.volume() / 2.0)
        B0, B0prime, E0, V0 = murnaghan(volumes, energies)
        relaxed_energy = E0 * 2.0
        # relax by minimizing stress, keeping box coordinates of ions fixed
        ase_tools.minimize_stress(system, 'BFGSLineSearch', 1e-4)
        self.assertAlmostEqual(relaxed_energy, system.energy(), places=5)
        # peturb ions and box, then restore by minimizing forces and stress 
        tm = np.array([[ 0.90, -0.03,  0.05],
                       [-0.03,  0.98,  0.04],
                       [ 0.05,  0.04,  1.07]])
        system.move_ions(tm.dot(np.array(system.ions_xyz_coords('a')).T).T, 'a')
        system.set_box(tm.dot(np.array(system.box_vectors('a')).T).T, 'a')
        ase_tools.minimize_forces_stress(system, 'BFGSLineSearch', 1e-4)
        self.assertAlmostEqual(relaxed_energy, system.energy(), places=5)

if __name__ == '__main__':
    unittest.main()
