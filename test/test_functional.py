# Copyright (c) 2019-2020 William C. Witt
# SPDX-License-Identifier: MIT

import math
import numpy as np
import os
import scipy.special
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), '../build/external/deft/'))
import pydeft as deft
sys.path.append(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), '../build/'))
import extras as profess
import tools_for_tests as tools

class TestFunctional(unittest.TestCase):

    def test_hartree(self):

        # set grid shape and box
        shape = (65, 68, 71)
        box = deft.Box([[25.0,0.0,0.0], [0.0,26.0,0.0], [0.0,0.0,27.0]])
        # construct test density (difference of two gaussians)
        # integrates to zero by design; if not, the hartree
        # potential explodes for periodic boundary conditions
        density = deft.Double3D(shape)
        w1 = 2.0
        w2 = 1.0
        def den(x,y,z):
            return (w1/np.pi)**1.5*np.exp(-w1*(x*x+y*y+z*z)) \
                    - (w2/np.pi)**1.5*np.exp(-w2*(x*x+y*y+z*z))
        density[...] = tools.get_function_on_grid(den, shape, box.vectors())
        # get exact hartree energy and exact hartree potential
        exact_ene = 1.0/np.sqrt(np.pi)*(
                np.sqrt(w1/2.0)-2.0*np.sqrt(w1*w2/(w1+w2))+np.sqrt(w2/2.0))
        def pot(x,y,z):
            r = np.sqrt(x*x+y*y+z*z)
            mask = r>1e-5
            pot = np.empty(r.shape)
            pot[np.invert(mask)] = 2.0*np.sqrt(w1/np.pi)-2.0*np.sqrt(w2/np.pi)
            pot[mask] = (scipy.special.erf(np.sqrt(w1)*r[mask])/r[mask]
                            - scipy.special.erf(np.sqrt(w2)*r[mask])/r[mask])
            return pot
        exact_pot = tools.get_function_on_grid(pot, shape, box.vectors())
        # compute energy and potential on grid; compare with exact quantities
        energy, potential = profess.Hartree(box).energy_potential(density)
        self.assertAlmostEqual(energy, exact_ene, places=6)
        self.assertTrue(np.allclose(potential[...], exact_pot, atol=1e-4))
        # set grid and box for testing stress
        # the grid is extremely coarse by design: there are subtleties
        # with the stress calculation arising from the half-complex 
        # storage format of the fft. if these subtleties aren't handled
        # correctly, this example fails by design.
        shape = (2, 2, 2)
        box = deft.Box([[ 3.0, 1.0, 1.0],
                        [-1.0, 4.0, 0.0],
                        [-0.5, 0.5, 4.0]])
        density = deft.Double3D(shape)
        density[...] = tools.get_function_on_grid(den, shape, box.vectors())
        stress = profess.Hartree(box).stress(density)
        def get_energy(box_vectors):
            new_box = deft.Box(box_vectors)
            return profess.Hartree(new_box).energy(
                    density*box.volume()/new_box.volume())
        approx_stress = tools.approximate_stress(get_energy, box.vectors())
        self.assertTrue(np.allclose(stress, approx_stress, atol=1e-10))

    def test_ion_electron(self):

        # choose grid and box
        shape = (49, 25, 25)
        box = deft.Box([[10, 0, 0],
                        [ 0, 5, 0],
                        [ 0, 0, 5]])
        # create test density
        density = deft.Double3D(shape)
        def harmonic_oscillator_3d_s_density(w,n,x,y,z):
            rr = x*x + y*y + z*z
            laguerre = scipy.special.eval_genlaguerre(n, 0.5, w*rr)
            return ((w/np.pi)**1.5
                        * (2**n*scipy.special.factorial(n)
                               / scipy.special.factorial2(2*n+1))
                        * laguerre * laguerre * np.exp(-w*rr))
        density[...] = tools.get_function_on_grid(
            lambda x,y,z: harmonic_oscillator_3d_s_density(7,0,x,y,z) \
                            + harmonic_oscillator_3d_s_density(7,1,x,y,z),
            shape,
            box.vectors(),
            [2.5, 0, 0])
        density[...] += tools.get_function_on_grid(
            lambda x,y,z: harmonic_oscillator_3d_s_density(5,0,x,y,z) \
                            + harmonic_oscillator_3d_s_density(5,1,x,y,z),
            shape,
            box.vectors(),
            [7.5, 0, 0])
        # compute and test energy
        ions = profess.Ions()
        r1 = 2.5
        r2 = 3.5
        ions.add_ion_type_harmonic_compactified(7.0, r1, r2)
        ions.add_ion_type_harmonic_compactified(5.0, r1, r2)
        xyz_coords = np.array([[2.5, 0.0, 0.0],
                               [7.5, 0.0, 0.0]])
        ion_type_ids = np.array([0, 1], dtype='int')
        ions.set_ions(xyz_coords, ion_type_ids)
        exact_energy = 30.0
        energy, potential = \
            profess.IonElectron(shape, box, ions).energy_potential(density)
        self.assertAlmostEqual(energy, exact_energy, places=2)
        energy, potential = \
            profess.IonElectron(shape, box, ions, 12).energy_potential(density)
        self.assertAlmostEqual(energy, exact_energy, places=2)
        energy, potential = \
            profess.IonElectron(shape, box, ions, 16).energy_potential(density)
        self.assertAlmostEqual(energy, exact_energy, places=2)
        # compute and test stress
        vectors = np.asarray(box.vectors())
        box_coords = np.linalg.inv(vectors.T).dot(xyz_coords.T).T
        vectors[1,0] = vectors[2,0] = vectors[2,1] = 1.0
        xyz_coords = (vectors.T.dot(box_coords.T)).T
        density *= (box.volume() / np.linalg.det(vectors))
        box.set(vectors)
        stress = profess.IonElectron(shape, box, ions).stress(density)
        def get_energy(box_vectors):
            new_box = deft.Box(box_vectors)
            new_loc = (box_vectors.T.dot(box_coords.T)).T
            ions.set_ions(box_vectors.T.dot(box_coords.T).T, ion_type_ids)
            f = profess.IonElectron(shape, new_box, ions)
            ions.set_ions(xyz_coords, ion_type_ids) # restore to unpeturbed
            return f.energy(density*box.volume()/new_box.volume())
        approx_stress = tools.approximate_stress(get_energy, box.vectors())
        self.assertTrue(np.allclose(stress, approx_stress, atol=1e-10))
        # compute and test forces
        forces = profess.IonElectron(shape, box, ions).forces(density)
        def get_energy(new_xyz_coords):
            ions.set_ions(new_xyz_coords, ion_type_ids)
            f = profess.IonElectron(shape, box, ions)
            ions.set_ions(xyz_coords, ion_type_ids) # restore to unpeturbed
            return f.energy(density)
        approx_forces = \
            tools.approximate_forces(get_energy, xyz_coords, box.vectors())
        self.assertTrue(np.allclose(forces, approx_forces, atol=0, rtol=1e-4))

    def test_libxc(self):

        # choose grid and box
        g = (19,19,19)
        box = deft.Box(8.0* np.eye(3))
        density = deft.Double3D(g)
        # create gaussian density
        w = 2.0
        def gaussian(x,y,z):
            return (w/np.pi)**1.5*np.exp(-w*(x*x+y*y+z*z))
        density[...] = tools.get_function_on_grid(gaussian, g, box.vectors())
        box = deft.Box(box.vectors())
        # ---------- LDA ----------
        # compute energy/potential
        energy, potential = \
            profess.PerdewZunger(box).energy_potential(density)
        # test energy/potential
        #exact_energy 
        #self.assertAlmostEqual(energy, exact_energy, places=6)
        #exact_potential = 
        #self.assertTrue(np.allclose(potential.data[...], exact_potential))
        # compute stress after peturbing box
        new_vectors = np.asarray(box.vectors())
        new_vectors[1:2,0] = 2.0
        new_vectors[0,1] = 3.0
        new_vectors[1,2] = 5.0
        density *= box.volume()/np.linalg.det(new_vectors)
        box.set(new_vectors)
        stress = profess.PerdewZunger(box).stress(density)
        # test stress
        def get_energy(box_vectors):
            new_box = deft.Box(box_vectors)
            return profess.PerdewZunger(new_box).energy(
                density*box.volume()/new_box.volume())
        approx_stress = tools.approximate_stress(get_energy, box.vectors())
        self.assertTrue(np.allclose(stress, approx_stress, atol=1e-10))
        # ---------- PBE ----------
        # compute energy/potential
        energy, potential = \
            profess.PerdewBurkeErnzerhof(box).energy_potential(density)
        # test energy/potential
        #exact_energy 
        #self.assertAlmostEqual(energy, exact_energy, places=6)
        #exact_potential = 
        #self.assertTrue(np.allclose(potential.data[...], exact_potential))
        # compute stress after peturbing box
        new_vectors = np.asarray(box.vectors())
        new_vectors[1:2,0] = 2.0
        new_vectors[0,1] = 3.0
        new_vectors[1,2] = 5.0
        density *= box.volume()/np.linalg.det(new_vectors)
        box.set(new_vectors)
        stress = profess.PerdewBurkeErnzerhof(box).stress(density)
        # test stress
        def get_energy(box_vectors):
            new_box = deft.Box(box_vectors)
            return profess.PerdewBurkeErnzerhof(new_box).energy(
                density*box.volume()/new_box.volume())
        approx_stress = tools.approximate_stress(get_energy, box.vectors())
        self.assertTrue(np.allclose(stress, approx_stress, atol=1e-10))

    def test_kinetic_class_a(self):

        # construct test density
        shape = (5, 5, 5)
        box = deft.Box([[1.0, 0.0, 0.0],
                        [0.7, 0.7, 0.0],
                        [0.8, 0.8, 0.8]])
        w = 3.0
        density = deft.array_from_lattice_sum(
            shape,
            box,
            np.array([[0,0,0]]),
            lambda x,y,z: np.exp(-(x*x+y*y+z*z)/(4.0*w)))
        # choose functional parameters
        a = 1.0
        b = 1.0
        def f(x): return np.exp(x) # 1.0 + x
        def fp(x): return np.exp(x) # 1.0
        # test potential
        den0 = deft.integrate(density,box) / box.volume()
        energy, potential = profess.KineticClassA(
            box,shape,a,b,f,fp,den0,False).energy_potential(density)
        def get_energy(density):
            return profess.KineticClassA(
                box,shape,a,b,f,fp,den0,False).energy(density)
        approx_pot = \
            tools.approximate_potential(get_energy, density, box.volume())
        self.assertTrue(
            np.allclose(potential[...], approx_pot, atol=0, rtol=1e-4))
        # test stress, case when den0 is held fixed as box changes
        den0 = 1.0
        stress = profess.KineticClassA(
            box,shape,a,b,f,fp,den0,True).stress(density)
        def get_energy(box_vectors):
            new_box = deft.Box(box_vectors)
            return profess.KineticClassA(
                new_box,shape,a,b,f,fp,den0,True
                ).energy(density*box.volume()/new_box.volume())
        approx_stress = tools.approximate_stress(get_energy, box.vectors())
        self.assertTrue(np.allclose(stress, approx_stress, atol=0, rtol=1e-7))
        # test stress, case when den0 is scaled with changing box
        den0 = deft.integrate(density,box) / box.volume()
        stress = profess.KineticClassA(
            box,shape,a,b,f,fp,den0,False).stress(density)
        def get_energy(box_vectors):
            new_box = deft.Box(box_vectors)
            den0 = deft.integrate(density,box) / new_box.volume()
            return profess.KineticClassA(
                new_box,shape,a,b,f,fp,den0,False
                ).energy(density*box.volume()/new_box.volume())
        approx_stress = tools.approximate_stress(get_energy, box.vectors())
        self.assertTrue(np.allclose(stress, approx_stress, atol=0, rtol=1e-7))

    def test_thomas_fermi(self):

        # choose grid and box
        shape = (12, 12, 12)
        box = deft.Box(5 * np.eye(3))
        density = deft.Double3D(shape)
        # create gaussian density and compute energy/potential
        w = 2.0
        def gaussian(x,y,z):
            return (w/np.pi)**1.5*np.exp(-w*(x*x+y*y+z*z))
        density[...] = \
            tools.get_function_on_grid(gaussian, shape, box.vectors())
        energy, potential = profess.ThomasFermi(box).energy_potential(density)
        # test energy/potential
        c0 = 3.0/10.0*(3.0*np.pi*np.pi)**(2.0/3.0)
        exact_energy = (3.0/5.0)**(3.0/2.0)/np.pi*c0*w
        self.assertAlmostEqual(energy, exact_energy, places=6)
        exact_potential = 5.0/3.0*c0*density[...]**(2.0/3.0)
        self.assertTrue(np.allclose(potential[...], exact_potential))
        # compute stress
        stress = profess.ThomasFermi(box).stress(density)
        # test stress
        def get_energy(box_vectors):
            new_box = deft.Box(box_vectors)
            return profess.ThomasFermi(new_box).energy(
                    density*box.volume()/new_box.volume())
        approx_stress = tools.approximate_stress(get_energy, box.vectors())
        self.assertTrue(np.allclose(stress, approx_stress, atol=1e-10))

    def test_weizsaecker(self):

        # choose grid and box
        shape = (19, 19, 19)
        box = deft.Box(8.0*np.eye(3))
        # create gaussian density and compute energy/potential
        w = 2.0
        def f(x,y,z):
            return (w/np.pi)**1.5*np.exp(-w*(x*x+y*y+z*z))
        density = deft.Double3D(shape)
        density[...] = tools.get_function_on_grid(f, shape, box.vectors())
        energy, potential = profess.Weizsaecker(box).energy_potential(density)
        # test energy
        exact_energy = 3.0 * w / 4.0
        self.assertAlmostEqual(energy, exact_energy, places=6)
        # test potential, comparing pot*sqrt(den) to damp error in vacuum
        def f(x,y,z):
            return -0.5*w*(w*(x*x+y*y+z*z)-3.0)
        exact_potential = tools.get_function_on_grid(f, shape, box.vectors())
        self.assertTrue(np.allclose(
            potential[...]*np.sqrt(density[...]), 
            exact_potential*np.sqrt(density[...]), atol=1e-6))
        # compute stress after peturbing box
        vectors = np.asarray(box.vectors())
        vectors[0,1:2] = 2.0
        vectors[1,0] = 3.0
        vectors[2,1] = 5.0
        density *= box.volume() / np.linalg.det(vectors)
        box.set(vectors)
        stress = profess.Weizsaecker(box).stress(density)
        # test stress
        def get_energy(box_vectors):
            new_box = deft.Box(box_vectors)
            return profess.Weizsaecker(new_box).energy(
                density*box.volume()/new_box.volume())
        approx_stress = tools.approximate_stress(get_energy, box.vectors())
        self.assertTrue(np.allclose(stress, approx_stress, atol=1e-10))

    def test_thomas_fermi_weizsaecker(self):

        # functional parameters
        a = 1.0
        b = 1.0
        tiny = 1e-12
        # choose grid and box
        shape = (35,35,35)
        box = deft.Box(8.0*np.eye(3))
        # create gaussian density and compute energy/potential
        w = 2.0
        def f(x,y,z):
            return (w/np.pi)**1.5*np.exp(-w*(x*x+y*y+z*z))
        density = deft.Double3D(shape)
        density[...] = tools.get_function_on_grid(f, shape, box.vectors())
        energy, potential = profess.ThomasFermiWeizsaecker(
            box, a, b, tiny).energy_potential(density)
        # test energy
        the_energy = 2.3495238708751
        self.assertAlmostEqual(energy, the_energy, places=8)
        # test potential, comparing pot*sqrt(den) to damp error in vacuum
        def f(x,y,z):
            return (0.5*3.0**(2.0/3.0)*np.pi**(1.0/3.0)*w
                        *np.exp(-2.0/3.0*w*(x*x+y*y+z*z))
                    -0.5*w*(w*(x*x+y*y+z*z)-3.0))
        the_potential = tools.get_function_on_grid(f, shape, box.vectors())
        self.assertTrue(np.allclose(
            potential[...]*density[...],
            the_potential*density[...], atol=1e-8, rtol=0))
        # compute stress after peturbing box
        vectors = np.asarray(box.vectors())
        vectors[0,1:2] = 2.0
        vectors[1,0] = 3.0
        vectors[2,1] = 5.0
        density *= box.volume() / np.linalg.det(vectors)
        box.set(vectors)
        stress = profess.ThomasFermiWeizsaecker(box,a,b,tiny).stress(density)
        # test stress
        def get_energy(box_vectors):
            new_box = deft.Box(box_vectors)
            return profess.ThomasFermiWeizsaecker(new_box,a,b,tiny).energy(
                density*box.volume()/new_box.volume())
        the_stress = tools.approximate_stress(get_energy, box.vectors())
        self.assertTrue(np.allclose(stress, the_stress, atol=1e-10, rtol=0))

    def test_luo_karasiev_trickey(self):

        # choose grid and box
        shape = (35,35,35)
        box = deft.Box(8.0*np.eye(3))
        # create gaussian density and compute energy/potential
        w = 2.0
        def f(x,y,z):
            return (w/np.pi)**1.5*np.exp(-w*(x*x+y*y+z*z))
        density = deft.Double3D(shape)
        density[...] = tools.get_function_on_grid(f, shape, box.vectors())
        energy, potential = profess.LuoKarasievTrickey(
            box).energy_potential(density)
        # test energy
        the_energy = 2.07914244900
        self.assertAlmostEqual(energy, the_energy, places=8)
        # test potential, comparing pot*sqrt(den) to damp error in vacuum
        def f(x,y,z):
            return (0.5*3.0**(2.0/3.0)*np.pi**(1.0/3.0)*w
                        *np.exp(-2.0/3.0*w*(x*x+y*y+z*z))
                    -0.5*w*(w*(x*x+y*y+z*z)-3.0))
        the_potential = tools.get_function_on_grid(f, shape, box.vectors())
        #self.assertTrue(np.allclose(
        #    potential[...]*density[...],
        #    the_potential*density[...], atol=1e-8, rtol=0))
        # compute stress after peturbing box
        vectors = np.asarray(box.vectors())
        vectors[0,1:2] = 2.0
        vectors[1,0] = 3.0
        vectors[2,1] = 5.0
        density *= box.volume() / np.linalg.det(vectors)
        box.set(vectors)
        stress = profess.LuoKarasievTrickey(box).stress(density)
        # test stress
        def get_energy(box_vectors):
            new_box = deft.Box(box_vectors)
            return profess.LuoKarasievTrickey(new_box).energy(
                density*box.volume()/new_box.volume())
        the_stress = tools.approximate_stress(get_energy, box.vectors())
        self.assertTrue(np.allclose(stress, the_stress, atol=1e-10, rtol=0))

if __name__ == '__main__':
    unittest.main()
