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

        .def(
            py::init<std::array<size_t,3>>(),
            R"(
                init docstring
            )")

        .def(
            "energy",
            &System::energy,
            R"(
                energy(unit='h')
                
                Compute total energy.
                
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

                Compute total enthalpy.
                
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
                forces(unit=?)
                
                Compute forces.
                
                Parameters
                ----------
                unit : str
                    Force unit.
                
                Returns
                -------
                forces : array_like
            )")

        .def(
            "stress",
            &System::stress,
            R"(
                stress(unit='h/b3')

                Compute stress tensor.

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
                
                Compute total pressure.
                
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
                
                Box vectors (cell vectors).
                
                Parameters
                ----------
                unit : str
                    Length unit.
                
                Returns
                -------
                box_vectors : array_like.
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
                grid_shape : array_like
             )")

        .def(
            "energy_cutoff",
            &System::energy_cutoff,
            R"(
                energy_cutoff(unit='h')
                
                Plane wave cutoff for the current box vectors and grid.
                
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

                For the ion potential defined in `filename`, add ions to the
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
                    Cartesian coordinates of ions.
                unit : str
                    Length unit for coords.

                Returns
                -------
                system : profess.System
            )",
            py::arg("z"),
            py::arg("coords"),
            py::arg("units")=std::string{"b"})

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
                add_luo_karasiev_trickey_functional

                Add Luo-Karasiev-Trickey kinetic energy functional.

                Parameters
                ----------
                a : float
                tiny_den : float

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
                    libxc ids

                Returns
                -------
                system : profess.System
            )")

        .def(
            "add_kinetic_class_a_functional",
            &System::add_kinetic_class_a_functional,
            R"(
                add_kinetic_class_a_functional(a,b,f,fp,den0=-1)

                Add a generic kinetic energy functional of type A.

                Parameters
                ----------
                a : float
                    exponent a
                b : float
                    exponent b
                f : function
                    f(x)
                fp : function
                    f'(x)
                den0 : float, optional
                    uniform reference density

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

                Add Perdew-Zunger local density approximation exchange-correlation functional.

                Parameters
                ----------

                Returns
                -------
            )")

        .def(
            "add_perrot_functional",
            &System::add_perrot_functional,
            R"(
                add_perrot_functional

                Add Perrot kinetic energy functional.

                Parameters
                ----------
                den0 : float
                    uniform reference density

                Returns
                -------
                system : profess.System
            )",
            py::arg("den0")=-1)

        .def(
            "add_smargiassi_madden_functional",
            &System::add_smargiassi_madden_functional,
            R"(
                add_smargiassi_madden_functional

                Add Smargiassi-Madden kinetic energy functional.

                Parameters
                ----------
                den0 : float
                    uniform reference density

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
                add_wang_teter_functional

                Add Wang-Teter kinetic energy functional.

                Parameters
                ----------
                den0 : float
                    uniform reference density

                Returns
                -------
            )",
            py::arg("den0")=-1)

        .def(
            "add_wang_govind_carter_functional",
            &System::add_wang_govind_carter_functional,
            R"(
                add_wang_govind_carter_functional
                Parameters
                ----------

                Returns
                -------
            )",
            py::arg("den0"),
            py::arg("alpha")=(5.0+std::sqrt(5.0))/6.0,
            py::arg("beta")=(5.0-std::sqrt(5.0))/6.0,
            py::arg("gamma")=2.7)

        .def(
            "add_wang_govind_carter_1999_i_functional",
            &System::add_wang_govind_carter_1999_i_functional,
            R"(
                add_wang_govind_carter_1999_i_functional
                Parameters
                ----------

                Returns
                -------
            )",
            py::arg("den0")=-1)

        .def(
            "add_weizsaecker_functional",
            &System::add_weizsaecker_functional,
            R"(
                add_weizsaecker_functional
                Parameters
                ----------

                Returns
                -------
            )")

        .def(
            "remove_functional",
            &System::remove_functional)

        .def(
            "add_electrons",
            &System::add_electrons)

        .def(
            "set_box",
            &System::set_box,
             R"(
                set_box(vectors, unit)

                Set the simulation box
                
                    vectors (array): box vectors
                    unit (str): length unit
              )",
             py::arg("vectors"),
             py::arg("unit")=std::string{"b"})

        .def(
            "move_ions",
            &System::move_ions,
            py::arg("xyz_coords"),
            py::arg("unit")=std::string{"b"})

        .def(
            "minimize_energy",
            &System::minimize_energy,
            "varies the electron density to minimize the system's energy",
            py::arg("energy_tol") = 1e-7,
            py::arg("window_size") = 3,
            py::arg("max_iter") = 1000)

        .def(
            "minimize_energy_tpsd", &System::minimize_energy_tpsd,
            "varies the electron density to minimize the system's energy",
            py::arg("energy_tol") = 1e-7,
            py::arg("window_size") = 3,
            py::arg("max_iter") = 1000)

        // -------------
        // class methods
        // -------------

        .def_static(
            "get_shape",
            &System::get_shape,
            "Docstring for get_shape");
}

}
