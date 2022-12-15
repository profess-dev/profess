// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "tools.hpp"
#include "minimizers.hpp"

#include "hartree.hpp"
#include "ion_electron.hpp"
#include "kinetic_class_a.hpp"
#include "luo_karasiev_trickey.hpp"
#include "perdew_burke_ernzerhof.hpp"
#include "perdew_zunger.hpp"
#include "perrot.hpp"
#include "smargiassi_madden.hpp"
#include "thomas_fermi.hpp"
#include "thomas_fermi_weizsaecker.hpp"
#include "wang_govind_carter_1999_i.hpp"
#include "wang_teter.hpp"
#include "weizsaecker.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace profess
{

namespace py = pybind11;

// helper function for instantiating a Functional of type F
template<typename F, typename ...ConstructorArgTypes>
void define_functional(py::module &m, std::string name)
{
    py::class_<F>(m, name.c_str())
        .def(py::init<ConstructorArgTypes...>())
        .def("energy", &F::energy)
        .def("energy_potential", &F::energy_potential)
        .def("stress", &F::stress);
}

PYBIND11_MODULE(extras, m) {

    m.def("interpolate_cubic_lagrange",
        py::overload_cast<std::vector<double>,double,double,double>(
                &interpolate_cubic_lagrange));
    m.def("interpolate_derivative_cubic_lagrange",
        py::overload_cast<std::vector<double>,double,double,double>(
            &interpolate_derivative_cubic_lagrange));

    m.def("two_point_steepest_descent",
          &two_point_steepest_descent,
          py::arg("function"),
          py::arg("x"),
          py::arg("f_tol")=1e-7,
          py::arg("f_window")=2,
          py::arg("g_tol")=1e-5,
          py::arg("max_iter")=100);

    py::class_<Ions>(m, "Ions")

        .def(py::init<>())

        .def("count", &Ions::count)
        .def("count_types", &Ions::count_types)
        .def("count_of_type", &Ions::count_of_type)
        .def("count_charge", &Ions::count_charge)

        .def("xyz_coords", py::overload_cast<>(&Ions::xyz_coords))
        .def("xyz_coords", py::overload_cast<>(&Ions::xyz_coords, py::const_))

        .def("type_ids", &Ions::type_ids)

        .def("xyz_coords_of_type", &Ions::xyz_coords_of_type)

        .def("add_ion_type_coulomb", &Ions::add_ion_type_coulomb)
        .def("add_ion_type_recpot", &Ions::add_ion_type_recpot)
        .def("add_ion_type_harmonic_compactified",
            &Ions::add_ion_type_harmonic_compactified)
        .def("add_ion_type_generic", &Ions::add_ion_type_generic)

        .def("set_ions", &Ions::set_ions);

    // Hartree
    define_functional<Hartree,Box>(m, "Hartree");
    // IonElectron
    py::class_<IonElectron>(m, "IonElectron")
        .def(
            py::init<std::array<size_t,3>, Box, Ions, int>(),
            py::arg("grid_shape"),
            py::arg("box"),
            py::arg("ions"),
            py::arg("spline_order")=-1)
        .def("energy", &IonElectron::energy)
        .def("energy_potential", &IonElectron::energy_potential)
        .def("stress", &IonElectron::stress)
        .def("forces", &IonElectron::forces);
    define_functional<KineticClassA,
        Box, std::array<size_t,3>,
        double, double,
        std::function<double(double)>, std::function<double(double)>,
        double, bool>(m, "KineticClassA");
    define_functional<LibXC,Box,std::vector<int>>(m, "LibXC");
    define_functional<LuoKarasievTrickey,Box>(m, "LuoKarasievTrickey");
    define_functional<PerdewBurkeErnzerhof,Box>(m, "PerdewBurkeErnzerhof");
    define_functional<PerdewZunger,Box>(m, "PerdewZunger");
    define_functional<ThomasFermi,Box>(m, "ThomasFermi");
    define_functional<ThomasFermiWeizsaecker,Box,double,double,double>(m, "ThomasFermiWeizsaecker");
    define_functional<Weizsaecker,Box>(m, "Weizsaecker");
    define_functional<Perrot,Box,std::array<size_t,3>,double,bool>(m, "Perrot");
//    define_functional<WangTeter,Box,double>(m, "WangTeter");
//    define_functional<SmargiassiMadden,Box,double>(m, "SmargiassiMadden");
//    define_functional<WangGovindCarter1999I,Box,double>(
//        m, "WangGovindCarter1999I");

}

}
