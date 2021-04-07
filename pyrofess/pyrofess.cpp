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

               compute system energy.

               Args:
                   unit (str): energy unit

               Returns:
                   energy (float)
            )",
            py::arg("unit")=std::string{"h"})

        .def(
            "enthalpy",
            &System::enthalpy,
            R"(
               enthalpy(unit='h')

               compute system enthalpy.

               Args:
                   unit (str): energy unit

               Returns:
                   energy (float)
            )",
            py::arg("unit")=std::string{"h"})

        .def(
            "forces",
            &System::forces,
            R"(
               forces(unit=?)

               compute forces on ions.

               Args:
                   unit (?): ?

               Returns:
                   forces (2darray): 
            )")

        .def("stress",
             &System::stress,
             R"(
                stress(unit)

                compute stress tensor.

                Args:
                    unit (str)

                Returns:
                    stress (array)
             )",
             py::arg("unit")=std::string{"h/b3"})

        .def("pressure",
             &System::pressure,
             R"(
                pressure(unit)

                compute system pressure.

                Args:
                    unit (str)

                Returns:
                    pressure (float)
             )",
             py::arg("unit")=std::string{"h/b3"})

        .def("volume",
             &System::volume,
             R"(
                volume(unit)

                report current volume.

                Args:
                    unit (str)

                Returns:
                    volume (float)
             )",
             py::arg("unit")=std::string{"b3"})

        .def("box_vectors",
             &System::box_vectors,
             R"(
                box_vectors(unit='b')

                report current box vectors (cell vectors).

                Args:
                    unit (str)

                Returns:
                    box_vectors (numpy array)
             )",
             py::arg("unit")=std::string{"b"})

        .def_readonly(
            "grid_shape",
            &System::grid_shape,
            R"(
                 grid_shape()

                 report grid size along all three dimensions.

                 Returns:
                     grid_shape (numpy array)
             )")

        .def(
            "energy_cutoff",
            &System::energy_cutoff,
            R"(
               energy_cutoff(unit)

               report plane wave cutoff for the current grid.

               Returns:
                   energy_cutoff (float)
            )",
            py::arg("unit")=std::string{"h"})

        .def(
            "total_ion_charge",
            &System::total_ion_charge,
            R"(
                total ion charge doc.
            )")

        .def("ions_xyz_coords",
             &System::ions_xyz_coords,
             R"(
                ions_xyz_coords(unit)

                cartesian coordinates of ions.

                Args:
                    unit (str)

                Returns:
                    ions_xyz_coords (numpy array)
             )",
             py::arg("unit")=std::string{"b"})

        .def_readwrite(
            "electron_density",
            &System::electron_density,
            R"(
                electron_density()

                electron density on a grid.

                Args:
                    unit (str)

                Returns:
                    electron_density (numpy array)
            )")

        // -----------------
        // modify the system
        // -----------------

        .def("add_ions",
             &System::add_ions,
             py::arg("filename"),
             py::arg("coords"),
             py::arg("units")=std::string{"b"})

        .def("add_coulomb_ions",
             &System::add_coulomb_ions,
             py::arg("z"),
             py::arg("coords"),
             py::arg("units")=std::string{"b"})

        .def("add_hartree_functional", &System::add_hartree_functional)

        .def("add_huang_carter_functional",
             &System::add_huang_carter_functional,
             py::arg("den0"))

        .def("add_ion_electron_functional",
             &System::add_ion_electron_functional)

        .def("add_luo_karasiev_trickey_functional",
             &System::add_luo_karasiev_trickey_functional,
             py::arg("a")=1.3,
             py::arg("tiny_den")=1e-12)

        .def("add_libxc_functional",
             py::overload_cast<std::vector<int>>
                 (&System::add_libxc_functional))

        .def("add_kinetic_class_a_functional",
             &System::add_kinetic_class_a_functional,
             py::arg("a"),
             py::arg("b"),
             py::arg("f"),
             py::arg("fp"),
             py::arg("den0")=-1)

        .def("add_perdew_burke_ernzerhof_functional",
              &System::add_perdew_burke_ernzerhof_functional)

        .def("add_perdew_zunger_functional",
              &System::add_perdew_zunger_functional)

        .def("add_perrot_functional",
            &System::add_perrot_functional,
            py::arg("den0")=-1)

        .def("add_smargiassi_madden_functional",
             &System::add_smargiassi_madden_functional,
             py::arg("den0")=-1)

        .def("add_thomas_fermi_functional",
             &System::add_thomas_fermi_functional)

        .def("add_wang_teter_functional",
             &System::add_wang_teter_functional,
             py::arg("den0")=-1)

        .def("add_wang_govind_carter_functional",
             &System::add_wang_govind_carter_functional,
             py::arg("den0"),
             py::arg("alpha")=(5.0+std::sqrt(5.0))/6.0,
             py::arg("beta")=(5.0-std::sqrt(5.0))/6.0,
             py::arg("gamma")=2.7)

        .def("add_wang_govind_carter_1999_i_functional",
             &System::add_wang_govind_carter_1999_i_functional,
             py::arg("den0")=-1)

        .def("add_weizsaecker_functional",
             &System::add_weizsaecker_functional)

        .def("remove_functional",
             &System::remove_functional)

        .def("distribute_electrons_uniformly",
             &System::distribute_electrons_uniformly)

        .def("set_box",
             &System::set_box,
             R"(
                set_box(vectors, unit)

                Set the simulation box
                
                Args:
                    vectors (array): box vectors
                    unit (str): length unit
              )",
             py::arg("vectors"),
             py::arg("unit")=std::string{"b"})

        .def("move_ions",
             &System::move_ions,
             py::arg("xyz_coords"),
             py::arg("unit")=std::string{"b"})

        .def("minimize_energy", &System::minimize_energy,
             "varies the electron density to minimize the system's energy",
             py::arg("energy_tol") = 1e-7,
             py::arg("window_size") = 3,
             py::arg("max_iter") = 1000)

        .def("minimize_energy_tpsd", &System::minimize_energy_tpsd,
             "varies the electron density to minimize the system's energy",
             py::arg("energy_tol") = 1e-7,
             py::arg("window_size") = 3,
             py::arg("max_iter") = 1000)

        // -------------
        // class methods
        // -------------

        .def_static("get_shape",
                    &System::get_shape,
                    "Docstring for get_shape");
}

}
