// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "functional.hpp"
#include "ions.hpp"

#include <memory>

namespace profess
{

class System
{
public:

    System(std::array<size_t,3> grid_shape);

    static std::array<size_t,3> get_shape(
            std::array<std::array<double,3>,3> box_vectors,
            double energy_cutoff,
            std::array<std::string,2> units={"b","h"});

    // rule of five
    // copy construction/assignment disabled b/c class has unique ptrs
    ~System() = default;
    System(const System&) = delete;
    System& operator=(const System&) = delete;
    System(System&&) = default;
    System& operator=(System&&) = default;

    // basic characteristics
    std::array<size_t,3> grid_shape;
    Box box;
    Double3D electron_density;
    Ions ions;
    std::vector<std::unique_ptr<Functional>> functionals;

    std::array<std::array<double,3>,3> box_vectors(std::string unit);
    double box_volume(std::string unit);

    System& add_ions(
        std::string type,
        std::vector<std::array<double,3>> coords,
        std::string unit={"b"});
    System& add_coulomb_ions(
        double z,
        std::vector<std::array<double,3>> coords,
        std::string unit={"b"});
    System& add_harmonic_ions();

    std::vector<std::array<double,3>> ions_xyz_coords(std::string unit);

    System& add_hartree_functional();
    System& add_huang_carter_functional(double den0);
    System& add_ion_electron_functional();
    System& add_luo_karasiev_trickey_functional(double a=1.3,
                                                double tiny_den=1e-12);
    System& add_libxc_functional(std::vector<int>);
    System& add_kinetic_class_a_functional(
        double a,
        double b,
        std::function<double(double)> f,
        std::function<double(double)> fp,
        double den0=-1);
    System& add_perdew_burke_ernzerhof_functional();
    System& add_perdew_zunger_functional();
    System& add_perrot_functional(double den0=-1);
    System& add_smargiassi_madden_functional(double den0=-1);
    System& add_thomas_fermi_functional();
    System& add_wang_govind_carter_functional(
        double den0,
        double alpha=(5.0+std::sqrt(5.0))/6.0,
        double beta=(5.0-std::sqrt(5.0))/6.0,
        double gamma=2.7);
    System& add_wang_govind_carter_1999_i_functional(double den0=-1);
    System& add_wang_teter_functional(double den0=-1);
    System& add_weizsaecker_functional();

    System& remove_functional(std::string name);

    // basic information about the system
    double energy(std::string unit={"h"});
    std::tuple<double, Double3D> energy_potential(bool compute_ion_ion=true);
    std::vector<std::array<double,3>> forces();
    std::array<std::array<double,3>,3> stress(std::string unit={"h/b3"});
    double volume(std::string unit={"b3"});
    double pressure(std::string unit={"h/b3"});
    double enthalpy(std::string unit={"h"});
    double energy_cutoff(std::string unit={"h"});
    double total_ion_charge();

    // basic manipulations
    System& distribute_electrons_uniformly(const double electrons);
    System& move_ions(
            std::vector<std::array<double,3>> xyz_coords,
            std::string length_unit={"b"});
    System& set_box(
            std::array<std::array<double,3>,3> vectors,
            std::string length_unit={"b"});
    System& minimize_energy(
        double energy_tol=1e-7,
        size_t window_size=3,
        size_t max_iter=1000);
    System& minimize_energy_tpsd(
        double energy_tol=1e-7,
        size_t window_size=3,
        size_t max_iter=1000);
    System& minimize_forces(
        double energy_tol=1e-5,
        size_t window_size=3,
        size_t max_iter=100);
};

}
