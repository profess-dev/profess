// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "system.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace profess
{

namespace py = pybind11;

PYBIND11_MODULE(profess, m) {
    py::options options;
    options.disable_function_signatures();


    m.doc() = "profess module documentation";

    py::class_<System>(m, "System")

        .def_static(
            "create",
            &System::create,
            R"(
                create(box_vectors, energy_cutoff, units=['b','h'])
                
                Create empty system, using a plane wave energy cutoff to determine the grid dimensions.
                
                Parameters
                ----------
                box_vectors : array_like
                    Lattice vectors, provided in the rows of a 3x3 list or array.
                energy_cutoff : float
                    Maximum plane wave energy, sets grid resolution.
                units : list of str, optional
                    Units for previous arguments.

                Returns
                -------
                system : profess.System
            )",
            py::arg("box_vectors"),
            py::arg("plane_wave_cutoff"),
            py::arg("units")=std::array<std::string,2>({"b","h"}))

        .def_static(
            "create_from_grid_shape",
            &System::create_from_grid_shape,
            R"(
                create_from_grid_shape(box_vectors, grid_shape, unit='b')
                
                Create empty system with a specific grid dimensions.
                
                Parameters
                ----------
                box_vectors : array_like
                    Lattice vectors, provided in the rows of a 3x3 list or array.
                grid_shape : sequence of 3 ints
                    Grid dimensions along each lattice vector.
                unit : str, optional
                    Length unit for box_vectors.

                Returns
                -------
                system : profess.System
            )",
            py::arg("box_vectors"),
            py::arg("grid_shape"),
            py::arg("unit")=std::string{"b"})

        .def(
            "energy",
            &System::energy,
            R"(
                energy(unit='h')
                
                Total energy.
                
                Parameters
                ----------
                unit : str
                    Energy unit.

                Returns
                -------
                energy : float
            )",
            py::arg("unit")=std::string{"h"})

        .def(
            "enthalpy",
            &System::enthalpy,
            R"(
                enthalpy(unit='h')

                Total enthalpy.
                
                Parameters
                ----------
                unit : str
                    Enthalpy unit.
                
                Returns
                -------
                enthalpy : float
            )",
            py::arg("unit")=std::string{"h"})

        .def(
            "forces",
            &System::forces,
            R"(
                forces(unit="h/b")
                
                Forces on ions.
                
                Parameters
                ----------
                unit : str
                    Force unit.
                
                Returns
                -------
                forces : array_like
            )",
            py::arg("unit")=std::string{"h/b"})

        .def(
            "stress",
            &System::stress,
            R"(
                stress(unit='h/b3')

                Stress tensor.

                Parameters
                ----------
                unit : str
                    Stress unit.

                Returns
                -------
                stress : array_like
            )",
            py::arg("unit")=std::string{"h/b3"})

        .def(
            "pressure",
            &System::pressure,
            R"(
                pressure(unit='h/b3')
                
                Pressure.
                
                Parameters
                ----------
                unit : str
                    Stress unit.
                
                Returns
                -------
                pressure : float
            )",
            py::arg("unit")=std::string{"h/b3"})

        .def(
            "volume",
            &System::volume,
            R"(
                volume(unit='b3')
                
                Box volume.
                
                Parameters
                ----------
                unit : str
                    Volume unit.
                
                Returns
                -------
                volume : float
            )",
            py::arg("unit")=std::string{"b3"})

        .def(
            "box_vectors",
            &System::box_vectors,
            R"(
                box_vectors(unit='b')
                
                Box vectors (lattice vectors).
                
                Parameters
                ----------
                unit : str
                    Length unit.
                
                Returns
                -------
                box_vectors : array_like
                    Vectors as rows of 3x3 array.
            )",
            py::arg("unit")=std::string{"b"})

        .def_readonly(
            "grid_shape",
            &System::grid_shape,
            R"(
                grid_shape()
                
                Grid size along all three dimensions.
                
                Returns
                -------
                grid_shape : sequence of 3 ints
             )")

        .def(
            "energy_cutoff",
            &System::energy_cutoff,
            R"(
                energy_cutoff(unit='h')
                
                Plane wave energy cutoff for the current box vectors and grid.
                
                Parameters
                ----------
                unit : str
                    Energy unit.
                
                Returns
                -------
                energy_cutoff : float
            )",
            py::arg("unit")=std::string{"h"})

        .def(
            "total_ion_charge",
            &System::total_ion_charge,
            R"(
                total_ion_charge

                Total charge of all ions.

                Returns
                -------
                total_ion_charge : float
            )")

        .def(
            "ions_xyz_coords",
            &System::ions_xyz_coords,
            R"(
                ions_xyz_coords(unit='b')
                
                Cartesian coordinates of ions.
                
                Parameters
                ----------
                unit : str
                     Length unit.
                
                Returns
                -------
                ions_xyz_coords : array_like
                    Ion coordinates in shape Nx3, where N is the number of
                    ions.
            )",
            py::arg("unit")=std::string{"b"})

        .def_readwrite(
            "electron_density",
            &System::electron_density,
            R"(
                electron_density
                
                Electron density on the grid.
                
                Parameters
                ----------
                unit : str
                    Density unit.
                
                Returns
                -------
                electron_density : deft_array
            )")

        // -----------------
        // modify the system
        // -----------------

        .def(
            "add_ions",
            &System::add_ions,
            R"(
                add_ions(filename, coords, unit='frac')

                Add ions to the system.

                Using the ion potential in `filename`, add ions to the
                locations indicated by `coords`.

                Parameters
                ----------
                filename : str
                    Potential file.
                coords : array_like
                    Ion coordinates in shape Nx3, where N is the number of
                    ions.
                unit : str, optional
                    Length unit for `coords`. Can be 'frac' (the default),
                    indicating fractional coordinates, or a valid length unit,
                    indicating Cartesian coordinates.

                Returns
                -------
                system : profess.System
            )",
            py::arg("filename"),
            py::arg("coords"),
            py::arg("units")=std::string{"b"})

        .def(
            "add_coulomb_ions",
            &System::add_coulomb_ions,
            R"(
                add_coulomb_ions(z, coords, unit)

                Add Coulomb ions to the system.

                Parameters
                ----------
                z : float
                    Ion charge.
                coords : array_like
                    Ion coordinates in shape Nx3, where N is the number of
                    ions.
                unit : str, optional
                    Length unit for `coords`. Can be 'frac' (the default),
                    indicating fractional coordinates, or a valid length unit,
                    indicating Cartesian coordinates.

                Returns
                -------
                system : profess.System
            )",
            py::arg("z"),
            py::arg("coords"),
            py::arg("units")=std::string{"b"},
            py::arg("cutoff")=-1.0)

//        .def(
//            "add_harmonic_ions",
//            &System::add_harmonic_ions,
//            R"(
//                add_harmonic_ions
//
//                Add harmonic ions to the system.
//
//                Returns
//                -------
//                system : profess.System
//            )")

        .def(
            "add_hartree_functional",
            &System::add_hartree_functional,
            R"(
                add_hartree_functional

                Add Hartree functional.

                Returns
                -------
                system : profess.System
            )")

        .def(
            "add_huang_carter_functional",
            &System::add_huang_carter_functional,
            R"(
                add_huang_carter_functional

                Add Huang Carter kinetic energy functional.

                Parameters
                ----------
                den0 : float
                    uniform reference density

                Returns
                -------
                system : profess.System
            )",
            py::arg("den0"))

        .def(
            "add_ion_electron_functional",
            &System::add_ion_electron_functional,
            R"(
                add_ion_electron_functional

                Add ion-electron interaction functional.

                Returns
                -------
                system : profess.System
            )")

        .def(
            "add_luo_karasiev_trickey_functional",
            &System::add_luo_karasiev_trickey_functional,
            R"(
                add_luo_karasiev_trickey_functional(a=1.3, tiny_den=1e-12)

                Add Luo-Karasiev-Trickey kinetic energy functional.

                Parameters
                ----------
                a : float, optional
                tiny_den : float, optional

                Returns
                -------
                system : profess.System
            )",
            py::arg("a")=1.3,
            py::arg("tiny_den")=1e-12)

        .def(
            "add_libxc_functional",
            py::overload_cast<std::vector<int>>(&System::add_libxc_functional),
            R"(
                add_libxc_functional

                Add libxc exchange-correlation functional.

                Parameters
                ----------
                xc_func_ids : list of int
                    Individual contributions from libxc.

                Returns
                -------
                system : profess.System
            )")

        .def(
            "add_generic_nonlocal_a_functional",
            &System::add_generic_nonlocal_a_functional,
            R"(
                add_kinetic_nonlocal_a_functional(a,b,f,fp,den0=-1)

                Add nonlocal kinetic energy functional of generic type A.

                Parameters
                ----------
                a : float
                    Exponent `a`.
                b : float
                    Exponent `b`.
                f : function
                    The function f(x).
                fp : function
                    The derivative f'(x).
                den0 : float, optional
                    Uniform reference density. By default (whenever`den0`<0),
                    uses the system-averaged electron density.

                Returns
                -------
                system : profess.System
            )",
            py::arg("a"),
            py::arg("b"),
            py::arg("f"),
            py::arg("fp"),
            py::arg("den0")=-1)

        .def(
            "add_perdew_burke_ernzerhof_functional",
            &System::add_perdew_burke_ernzerhof_functional,
            R"(
                add_perdew_burke_ernzerhof_functional

                Add Perdew-Burke-Ernzerhof generalized gradient approximation exchange-correlation functional.

                Returns
                -------
                system : profess.System
            )")

        .def(
            "add_perdew_zunger_functional",
            &System::add_perdew_zunger_functional,
            R"(
                add_perdew_zunger_functional

                Add local density approximation functional,
                Perdew-Zunger parameterization.

                Returns
                -------
                system : profess.System
            )")

        .def(
            "add_perrot_functional",
            &System::add_perrot_functional,
            R"(
                add_perrot_functional(den0=-1)

                Add Perrot kinetic energy functional.

                Parameters
                ----------
                den0 : float, optional
                    Uniform reference density. By default (whenever`den0`<0),
                    uses the system-averaged electron density.

                Returns
                -------
                system : profess.System
            )",
            py::arg("den0")=-1)

        .def(
            "add_smargiassi_madden_functional",
            &System::add_smargiassi_madden_functional,
            R"(
                add_smargiassi_madden_functional(den0=-1)

                Add Smargiassi-Madden kinetic energy functional.

                Parameters
                ----------
                den0 : float, optional
                    Uniform reference density. By default (whenever`den0`<0),
                    uses the system-averaged electron density.

                Returns
                -------
                system : profess.System
            )",
            py::arg("den0")=-1)

        .def(
            "add_thomas_fermi_functional",
            &System::add_thomas_fermi_functional,
            R"(
                add_thomas_fermi_functional

                Add Thomas-Fermi kinetic energy functional.

                Returns
                -------
                system : profess.System
            )")

        .def(
            "add_wang_teter_functional",
            &System::add_wang_teter_functional,
            R"(
                add_wang_teter_functional(den0=-1)

                Add Wang-Teter kinetic energy functional.

                Parameters
                ----------
                den0 : float, optional
                    Uniform reference density. By default (whenever`den0`<0),
                    uses the system-averaged electron density.

                Returns
                -------
                system : profess.System
            )",
            py::arg("den0")=-1)

        .def(
            "add_wang_govind_carter_functional",
            &System::add_wang_govind_carter_functional,
            R"(
                add_wang_govind_carter_functional(den0=-1, gamma=2.7)

                Add Wang-Govind-Carter kinetic energy functional (density-dependent
                kernel).

                Parameters
                ----------
                den0 : float, optional
                    Uniform reference density. By default (whenever`den0`<0),
                    uses the system-averaged electron density.
                gamma : float, optional

                Returns
                -------
                system : profess.System
            )",
            py::arg("den0")=-1.0,
            py::arg("alpha")=(5.0+std::sqrt(5.0))/6.0,
            py::arg("beta")=(5.0-std::sqrt(5.0))/6.0,
            py::arg("gamma")=2.7)

        .def(
            "add_wang_govind_carter_1999_i_functional",
            &System::add_wang_govind_carter_1999_i_functional,
            R"(
                add_wang_govind_carter_1999_i_functional

                Add Wang-Govind-Carter functional (density-independent kernel).

                Parameters
                ----------
                den0 : float, optional
                    Uniform reference density. By default (whenever`den0`<0),
                    uses the system-averaged electron density.

                Returns
                -------
                system : profess.System
            )",
            py::arg("den0")=-1)

        .def(
            "add_weizsaecker_functional",
            &System::add_weizsaecker_functional,
            R"(
                add_weizsaecker_functional

                Add Weizsaecker kinetic energy functional.

                Returns
                -------
                system : profess.System
            )")

        .def(
            "remove_functional",
            &System::remove_functional)

        .def(
            "add_ion_ion_interaction",
            &System::add_ion_ion_interaction,
            R"(
                add_ion_ion_interaction()

                Add ion-ion interaction energy.

                Returns
                -------
                system : profess.System
            )")

        .def(
            "add_electrons",
            &System::add_electrons,
             R"(
                add_electrons(num=-1)

                Add electrons to the system.

                The electrons are distributed uniformly.
                
                Parameters
                ----------
                num : float, optional
                    By default (whenever `num` < 0), add a number of
                    electrons equal to the total ion charge in the
                    system.

                Returns
                -------
                system : profess.System
              )",
            py::arg("electrons")=-1.0)

        .def(
            "set_box",
            &System::set_box,
             R"(
                set_box(vectors, unit='h')

                Set the simulation box size and shape by changing the box vectors.

                The ion positions are scaled to preserve their fractional
                coordinates.

                Parameters
                ----------
                vectors : array_like
                    New box vectors, in rows of 3x3 list or array.
                unit : str, optional
                    Length unit for `vectors`.

                Returns
                -------
                system : profess.System
              )",
             py::arg("vectors"),
             py::arg("unit")=std::string{"b"})

        .def(
            "move_ions",
            &System::move_ions,
             R"(
                move_ions(xyz_coords, unit='b')

                Move the ions to new coordinates, keeping the box fixed.

                Parameters
                ----------
                xyz_coords : array_like
                    Ion coordinates in shape Nx3, where N is the number of
                    ions.
                unit : str, optional
                    Length unit.

                Returns
                -------
                system : profess.System
              )",
            py::arg("xyz_coords"),
            py::arg("unit")=std::string{"b"})

        .def(
            "minimize_energy",
            &System::minimize_energy,
             R"(
                minimize_energy(energy_tol=1e-7, window_size=3, max_iter=1000)

                Minimize total energy with LBFGS algorithm.

                Parameters
                ----------
                energy_tol : float, optional
                    Convergence tolerance.
                window_size : int, optional
                    Convergence achieved after the energy changes by less than
                    `energy_tol` over `window_size` steps.
                max_iter : int, optional
                    Maximum iterations permitted.

                Returns
                -------
                system : profess.System
              )",
            py::arg("energy_tol") = 1e-7,
            py::arg("window_size") = 3,
            py::arg("max_iter") = 1000)

        .def(
            "minimize_energy_tpsd",
            &System::minimize_energy_tpsd,
             R"(
                minimize_energy(energy_tol=1e-7, window_size=3, max_iter=1000)

                Minimize total energy with two-point steepest descent algorithm.

                Parameters
                ----------
                energy_tol : float, optional
                    Convergence tolerance.
                window_size : int, optional
                    Convergence achieved after the energy changes by less than
                    `energy_tol` over `window_size` steps.
                max_iter : int, optional
                    Maximum iterations permitted.

                Returns
                -------
                system : profess.System
              )",
            py::arg("energy_tol") = 1e-7,
            py::arg("window_size") = 3,
            py::arg("max_iter") = 1000);

}

}
