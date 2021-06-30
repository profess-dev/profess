# Copyright (c) 2019-2020 William C. Witt
# SPDX-License-Identifier: MIT

import numpy as np
import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), '../build/external/deft/'))
import pydeft as deft
sys.path.append(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), '../build/'))
import profess

class TestSystem(unittest.TestCase):

    def test_hydrogen_atom(self):

        system = (
            profess.System.create(20*np.identity(3), 10) 
            .add_coulomb_ions(1.0, [[0,0,0]], cutoff=15)
            .add_electrons()
            .add_weizsaecker_functional()
            .add_ion_electron_functional()
        )
        system.minimize_energy()
        self.assertAlmostEqual(system.energy(), -0.5, places=2)

    def test_fcc_aluminum(self):

        shape = (18,18,18)
        box_vectors = 4.050 * np.array([[0.5, 0.5, 0.0],
                                        [0.0, 0.5, 0.5],
                                        [0.5, 0.0, 0.5]])
        # w + tf
        system = (
            profess.System.create_from_grid_shape(box_vectors, shape, 'a')
            .add_ions('potentials/al.gga.recpot', np.zeros([1,3]), 'a')
            .add_electrons()
            .add_weizsaecker_functional()
            .add_thomas_fermi_functional()
            .add_hartree_functional()
            .add_ion_electron_functional()
            .add_perdew_burke_ernzerhof_functional()
            .add_ion_ion_interaction()
        )
        system.minimize_energy()
        self.assertAlmostEqual(system.energy('h'), -2.083724968, places=5)

        # wang-teter
        system = (
            profess.System.create_from_grid_shape(box_vectors, shape, 'a')
            .add_ions('potentials/al.gga.recpot', np.zeros([1,3]), 'a')
            .add_electrons()
            .add_wang_teter_functional()
            .add_hartree_functional()
            .add_perdew_burke_ernzerhof_functional()
            .add_ion_electron_functional()
            .add_ion_ion_interaction()
        )
        system.minimize_energy()

        # wang-govind-carter
        system = (
            profess.System.create_from_grid_shape(box_vectors, shape, 'a')
            .add_ions('potentials/al.gga.recpot', np.zeros([1,3]), 'a')
            .add_electrons()
            .add_wang_govind_carter_functional()
            .add_hartree_functional()
            .add_perdew_burke_ernzerhof_functional()
            .add_ion_electron_functional()
            .add_ion_ion_interaction()
        )
        system.minimize_energy()

    def test_bcc_lithium(self):
 
        box_vecs = 3.48 * np.eye(3)
        shape = (18,18,18)
        system = (
            profess.System.create_from_grid_shape(box_vecs, shape, 'a')
            .add_ions(
                'potentials/li.gga.recpot',
                box_vecs[0,0]*np.array([[0.0,0.0,0.0],[0.5,0.5,0.5]]),
                'a')
            .add_electrons()
            .add_smargiassi_madden_functional()
            .add_hartree_functional()
            .add_perdew_burke_ernzerhof_functional()
            .add_ion_electron_functional()
            .add_ion_ion_interaction()
        )
        system.minimize_energy()

if __name__ == '__main__':
    unittest.main()
