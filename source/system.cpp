// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "system.hpp"

#include "hartree.hpp"
#include "huang_carter.hpp"
#include "ion_electron.hpp"
#include "ion_ion.hpp"
#include "kinetic_class_a.hpp"
#include "luo_karasiev_trickey.hpp"
#include "libxc.hpp"
#include "perdew_burke_ernzerhof.hpp"
#include "perdew_zunger.hpp"
#include "perrot.hpp"
#include "smargiassi_madden.hpp"
#include "thomas_fermi.hpp"
#include "wang_govind_carter.hpp"
#include "wang_govind_carter_1999_i.hpp"
#include "wang_teter.hpp"
#include "weizsaecker.hpp"

#include "fourier.hpp"

#include "c_lbfgs.h"
#include "c_real_space_electrostatic_sum.h"

#include "minimizers.hpp"
#include "units.hpp"

#include <algorithm>
#include <limits>
#include <memory>

#include <iostream>

namespace profess {

System::System(std::array<size_t,3> grid_shape)
    : box({{{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}}}),
      grid_shape(grid_shape),
      electron_density(Double3D(grid_shape))
{}

System System::create(
    std::array<std::array<double,3>,3> box_vectors,
    double energy_cutoff,
    std::array<std::string,2> units)
{
    box_vectors = units::convert_length(box_vectors, units[0], {"b"});
    energy_cutoff = units::convert_energy(energy_cutoff, units[1], {"h"}); 
    double wavevector_cutoff = std::sqrt(2.0*energy_cutoff);
    auto a = box_vectors;
    double l0 = std::sqrt(a[0][0]*a[0][0]+a[0][1]*a[0][1]+a[0][2]*a[0][2]);
    double l1 = std::sqrt(a[1][0]*a[1][0]+a[1][1]*a[1][1]+a[1][2]*a[1][2]);
    double l2 = std::sqrt(a[2][0]*a[2][0]+a[2][1]*a[2][1]+a[2][2]*a[2][2]);
    size_t n0 = std::ceil(wavevector_cutoff*l0/(2.0*M_PI));
    size_t n1 = std::ceil(wavevector_cutoff*l1/(2.0*M_PI));
    size_t n2 = std::ceil(wavevector_cutoff*l2/(2.0*M_PI));
    System system({2*n0+1, 2*n1+1, 2*n2+1});
    system.set_box(box_vectors);
    return system;
}

System System::create_from_grid_shape(
    std::array<std::array<double,3>,3> box_vectors,
    std::array<size_t,3> grid_shape,
    std::string unit)
{
    System system(grid_shape);
    system.set_box(box_vectors, unit);
    return system;
}

std::array<std::array<double,3>,3> System::box_vectors(std::string unit)
{
    return units::convert_length(box.vectors(), {"b"}, unit);
}

double System::box_volume(std::string unit)
{
    return units::convert_volume(box.volume(), {"b3"}, unit);
}

System& System::add_ions(
    std::string filename,
    std::vector<std::array<double,3>> coords,
    std::string unit)
{
    ions.add_ion_type_recpot(filename);
    coords = units::convert_length(coords, unit, {"b"});
    auto ion_coords = ions.xyz_coords();
    ion_coords.insert(ion_coords.end(), coords.begin(), coords.end());
    auto ion_ids = ions.type_ids();
    size_t id = (ion_ids.size() == 0) ? 0 :
        *std::max_element(ion_ids.begin(), ion_ids.end()) + 1;
    auto new_ids = std::vector<size_t>(coords.size(), id);
    ion_ids.insert(ion_ids.end(), new_ids.begin(), new_ids.end());
    ions.set_ions(ion_coords, ion_ids); 
    return *this;
}

System& System::add_coulomb_ions(
    double z,
    std::vector<std::array<double,3>> coords,
    std::string unit,
    double cutoff)
{
    ions.add_ion_type_coulomb(z, cutoff);
    coords = units::convert_length(coords, unit, {"b"});
    auto ion_coords = ions.xyz_coords();
    ion_coords.insert(ion_coords.end(), coords.begin(), coords.end());
    auto ion_ids = ions.type_ids();
    size_t id = (ion_ids.size() == 0) ? 0 :
        *std::max_element(ion_ids.begin(), ion_ids.end()) + 1;
    auto new_ids = std::vector<size_t>(coords.size(), id);
    ion_ids.insert(ion_ids.end(), new_ids.begin(), new_ids.end());
    ions.set_ions(ion_coords, ion_ids); 
    return *this;
}

System& System::add_harmonic_ions()
{
    return *this;
}

std::vector<std::array<double,3>> System::ions_xyz_coords(std::string unit)
{
    return units::convert_length(ions.xyz_coords(), {"b"}, unit);
}

System& System::add_hartree_functional()
{
    functionals.emplace_back(std::make_unique<Hartree>(box));
    return *this;
}

System& System::add_huang_carter_functional(double den0)
{
    functionals.emplace_back(
        std::make_unique<HuangCarter>(box,grid_shape,den0));
    return *this;
}

System& System::add_ion_electron_functional()
{
    functionals.emplace_back(
            std::make_unique<IonElectron>(
                    electron_density.shape(),
                    box.vectors(),
                    ions));
    return *this;
}

System& System::add_luo_karasiev_trickey_functional(double a, double tiny_den)
{
    functionals.emplace_back(
        std::make_unique<LuoKarasievTrickey>(box, a, tiny_den));
    return *this;
}

System& System::add_libxc_functional(std::vector<int> ids)
{
    functionals.emplace_back(std::make_unique<LibXC>(box,ids));
    return *this;
}

System& System::add_generic_nonlocal_a_functional(
        double a,
        double b,
        std::function<double(double)> f,
        std::function<double(double)> fp,
        double den0)
{
    if (den0 < 0) den0 = total_ion_charge() / volume();
    functionals.emplace_back(std::make_unique<KineticClassA>(
            box, grid_shape, a, b, f, fp, den0));
    return *this;
}

System& System::add_perdew_burke_ernzerhof_functional() {
    functionals.emplace_back(std::make_unique<PerdewBurkeErnzerhof>(box));
    return *this;
}

System& System::add_perdew_zunger_functional() {
    functionals.emplace_back(std::make_unique<PerdewZunger>(box));
    return *this;
}

System& System::add_perrot_functional(double den0)
{
    if (den0 < 0) den0 = total_ion_charge() / volume();
    functionals.emplace_back(std::make_unique<Perrot>(
            box.vectors(), grid_shape, den0));
    return *this;
}

System& System::add_smargiassi_madden_functional(double den0)
{
    if (den0 < 0) den0 = total_ion_charge() / volume();
    functionals.emplace_back(std::make_unique<SmargiassiMadden>(
            box.vectors(), grid_shape, den0));
    return *this;
}

System& System::add_thomas_fermi_functional()
{
    functionals.emplace_back(std::make_unique<ThomasFermi>(box));
    return *this;
}

System& System::add_wang_govind_carter_functional(
    double den0, double alpha, double beta, double gamma)
{
    if (den0 < 0) den0 = total_ion_charge() / volume();
    functionals.emplace_back(
        std::make_unique<WangGovindCarter>(
            box.vectors(), grid_shape, den0, alpha, beta, gamma));
    return *this;
}

System& System::add_wang_govind_carter_1999_i_functional(double den0)
{
    if (den0 < 0) den0 = total_ion_charge() / volume();
    functionals.emplace_back(std::make_unique<WangGovindCarter1999I>(
            box.vectors(), grid_shape, den0));
    return *this;
}

System& System::add_wang_teter_functional(double den0)
{
    if (den0 < 0) den0 = total_ion_charge() / volume();
    functionals.emplace_back(std::make_unique<WangTeter>(
            box.vectors(), grid_shape, den0));
    return *this;
}

System& System::add_weizsaecker_functional()
{
    functionals.emplace_back(std::make_unique<Weizsaecker>(box));
    return *this;
}

System& System::remove_functional(std::string name)
{
    for (size_t i=0; i<functionals.size(); ++i) {

        if (functionals[i]->name() == name) functionals.erase(functionals.begin()+i);
    }
    return *this;
}

System& System::add_ion_ion_interaction()
{
    _ion_ion_interaction = true;
    return *this;
}

double System::energy(std::string unit)
{
    double energy;
    Double3D potential(electron_density);
    std::tie(energy, potential) = energy_potential();
    return units::convert_energy(energy, {"h"}, unit);
}

Double3D System::external_potential()
{
    std::cout << "WARNING: external_potential routine is still experimental." << std::endl;
    Double3D potential(electron_density);
    for (auto &f : functionals) {
        if (f->name() == "ion_electron") {
            double energy;
            std::tie(energy, potential) = f->energy_potential(electron_density);
            return potential;
        }
    }
    return potential.fill(0.0);
}

std::tuple<double, Double3D> System::energy_potential(bool compute_ion_ion)
{
    double energy, ene;
    Double3D potential(electron_density);
    Double3D pot(electron_density);

    energy = 0.0;
    potential.fill(0.0);

    for (auto &f : functionals) {
        std::tie(ene, pot) = f->energy_potential(electron_density);
        energy += ene;
        potential += pot;
        //std::cout << f->name() << "->  " << ene << std::endl;
    }

    if (_ion_ion_interaction) {
        if (compute_ion_ion) {
            energy += IonIon().energy(box.vectors(), ions);
        }
    }

    return std::make_tuple(energy, potential);
}

std::vector<std::array<double,3>> System::forces(std::string unit)
{
    profess::IonElectron IonElectron(
            electron_density.shape(),
            box.vectors(),
            ions);
    auto forces = IonElectron.forces(electron_density);
    const auto forces_ii = IonIon().forces(box.vectors(), ions);
    for (size_t i=0; i<forces.size(); ++i) {
        for (size_t j=0; j<3; ++j) {
            forces[i][j] += forces_ii[i][j];
        }
    }
    return units::convert_force(forces, {"h/b"}, unit);
}

std::array<std::array<double,3>,3> System::stress(std::string unit)
{
    std::array<std::array<double,3>,3> stress =
            {{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};

    for (auto &f : functionals) {
        auto s = f->stress(electron_density);
        for (size_t i=0; i<3; ++i) {
            for (size_t j=0; j<3; ++j) {
                stress[i][j] += s[i][j];
            }
        }
    }

    // ion-ion
    auto s = IonIon().stress(
                box.vectors(),
                ions);
    for (size_t i=0; i<3; ++i) {
        for (size_t j=0; j<3; ++j) {
            stress[i][j] += s[i][j];
        }
    }

    return units::convert_pressure(stress, {"h/b3"}, unit);
}

double System::volume(std::string unit)
{
    return units::convert_volume(box.volume(), {"b3"}, unit);
}

double System::pressure(std::string unit)
{
    auto s = stress();
    double p = -(s[0][0]+s[1][1]+s[2][2])/3.0;
    return units::convert_pressure(p, {"h/b3"}, unit);
}

double System::enthalpy(std::string unit)
{
    const double h = energy() + pressure()*volume();
    return units::convert_energy(h, {"h"}, unit);
}

double System::energy_cutoff(std::string unit)
{
    auto lengths = box.lengths();
    double e0 = 0.5*std::pow(2.0*M_PI*std::floor(grid_shape[0]/2)/lengths[0],2);
    double e1 = 0.5*std::pow(2.0*M_PI*std::floor(grid_shape[1]/2)/lengths[1],2);
    double e2 = 0.5*std::pow(2.0*M_PI*std::floor(grid_shape[2]/2)/lengths[2],2);
    double cutoff = std::min(e0, std::min(e1, e2));
    return units::convert_energy(cutoff, {"h"}, unit);
}

double System::total_ion_charge()
{
    return ions.count_charge();
}

System& System::add_electrons(double electrons)
{
    if (electrons < 0.0) {
        electron_density.fill(total_ion_charge() / box.volume());
    } else {
        electron_density.fill(electrons / box.volume());
    }
    return *this;
}

System& System::move_ions(
        std::vector<std::array<double,3>> xyz_coords,
        std::string length_unit)
{
    xyz_coords = units::convert_length(xyz_coords, length_unit, {"b"});

    // the following transforms the incoming xyz_coords to ensure
    // they lie within the box. in the future, may want this
    // to be an option, rather than performed by default.
    if (true) {

        // put box vectors in columns of matrix 'a'
        std::array<std::array<double,3>,3> a;
        for (size_t i=0; i<3; ++i) {
            for (size_t j=0; j<3; ++j) {
                a[i][j] = box.vectors()[j][i];
            }
        }

        // get inverse of box vector matrix
        std::array<std::array<double,3>,3> ainv;
        auto b = box.recip_vectors();
        for (size_t i=0; i<3; ++i) {
            for (size_t j=0; j<3; ++j) {
                // no transpose because b = [b1 b2 b3]^T
                ainv[i][j] = b[i][j] / (2.0*M_PI);
            }
        }

        // create function for matrix multiplication
        auto m_dot_v = [](std::array<std::array<double,3>,3> m, 
                          std::array<double,3> v) {
            std::array<double,3> r;
            r[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
            r[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
            r[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
            return r;
        };

        // ensure xyz_coords are within box
        for (size_t i=0; i<ions.count(); ++i) {
            auto coords = m_dot_v(ainv, xyz_coords[i]);
            for (size_t j=0; j<3; ++j) coords[j] -= std::floor(coords[j]);
            xyz_coords[i] = m_dot_v(a, coords);
        }
    }

    // save new coordinates
    ions.xyz_coords() = xyz_coords;

    // regenerate ion-electron functional
    for (auto &f : functionals) {
        if (f->name() == "ion_electron") {
            f = std::make_unique<IonElectron>(
                    grid_shape, box.vectors(), ions);
        }
    }

    return *this;    
}

System& System::set_box(
        std::array<std::array<double,3>,3> vectors,
        std::string length_unit)
{
    auto new_box = Box(units::convert_length(vectors, length_unit, {"b"}));

    // find inverse(A), assuming A=[a1 a2 a3]
    std::array<std::array<double,3>,3> ainv;
    auto b = box.recip_vectors();
    for (size_t i=0; i<3; ++i) {
        for (size_t j=0; j<3; ++j) {
            // no transpose here because b is stored as [b1 b2 b3]^T
            ainv[i][j] = b[i][j] / (2.0*M_PI);
        }
    }

    // create function for matrix multiplication
    auto m_dot_v = [](std::array<std::array<double,3>,3> m, 
                      std::array<double,3> v) {
        std::array<double,3> r;
        r[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
        r[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
        r[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
        return r;
    };
    
    // get box coordinates of ions
    const auto xyz_coords = ions.xyz_coords();
    std::vector<std::array<double,3>> box_coords(xyz_coords.size());
    for (size_t i=0; i<xyz_coords.size(); ++i) {
        box_coords[i] = m_dot_v(ainv, xyz_coords[i]);
    }

    // update system box
    electron_density *= box.volume()/new_box.volume();
    box = new_box;

    // put box vectors in columns of matrix 'a'
    std::array<std::array<double,3>,3> a;
    for (size_t i=0; i<3; ++i) {
        for (size_t j=0; j<3; ++j) {
            a[i][j] = box.vectors()[j][i];
        }
    }

    // update ion xyz coordinates
    for (size_t i=0; i<xyz_coords.size(); ++i) {
        ions.xyz_coords()[i] = m_dot_v(a, box_coords[i]);
    }

    // update functionals
    for (auto &f : functionals) f->set_box(box.vectors());

    // need a way for ion_ion to update the ions positions
    for (auto &f : functionals) {
        if (f->name() == "ion_electron") {
            f = std::make_unique<IonElectron>(
                    grid_shape, box.vectors(), ions);
        }
    }

    return *this;    
}

System& System::minimize_energy(
        double energy_tol,
        size_t window_size,
        size_t max_iter)
{
    // ----- define LBFGS variables, giving values to most -----
    int n = electron_density.size();       // problem size
    int m = 6;                             // iterations to store, 3<=m<=7
    Double3D x(electron_density);          // x array
    double f;                              // f(x)
    Double3D g(electron_density);          // df/dx array (gradient)
    int diagco = 0;                        // omit diagonal matrix Hk0
    std::vector<double> diag(n);           // diagonal matrix array
    std::array<int,2> iprint({-1,0});      // {1,0}->always print, but minimal
    double eps = 0.0;                      // stopping tolerance
    double xtol =                          // machine precision estimate
            std::numeric_limits<double>::epsilon();
    std::vector<double> w(n*(2*m+1)+2*m);  // work array
    int iflag = 0;                         // info flag

    double N = integrate(electron_density, box);
    x.compute_sqrt();

    std::vector<double> energy_history(window_size);
    std::iota(energy_history.begin(), energy_history.end(), 1.0);

    for (size_t i=0; i<max_iter; ++i) {

        // unpack density
        electron_density = x;
        electron_density *= electron_density;
        double Nbar = integrate(electron_density, box);
        electron_density *= (N/Nbar);

        // get functional potentials
        Double3D pot(electron_density);
        std::tie(f, pot) = energy_potential(false);
    
        // compute gradient
        g = pot;
        pot *= electron_density;
        double mu = 1.0 / N * integrate(pot, box);
        g -= mu;
        g *= x;
        g *= (2.0 * N / Nbar * box.volume() / electron_density.size());
    
        // conduct lbfgs step
        c_lbfgs_step(&n, &m, x.data(), &f, g.data(), &diagco, diag.data(), 
                     iprint.data(), &eps, &xtol, w.data(), &iflag);

        // check for lbfgs error
        if (iflag < 0) throw(std::runtime_error("lbfgs error"));

        // check for convergence
        if (std::all_of(
                energy_history.begin(),
                energy_history.end(),
                [energy_tol, f](double e){
                    return std::abs(e-f) < energy_tol;})) {
            break;
        } else {
            for (size_t j=window_size-1; j>0; --j) {
                energy_history[j] = energy_history[j-1];
            }
            energy_history[0] = f;
        }

        // check against max_iter
        if (i+1 == max_iter) throw(std::runtime_error("max_iter exceeded."));
    }
    return (*this);
}

System& System::minimize_energy_tpsd(
    double energy_tol,
    size_t window_size,
    size_t max_iter)
{
    double N = integrate(electron_density, box);
    std::vector<double> x(electron_density.size());
    for (size_t i=0; i<x.size(); ++i) x[i] = std::sqrt(electron_density(i));

    std::function<std::tuple<double,std::vector<double>>(
                      std::vector<double>)> compute_f_and_g =
        [this,N](std::vector<double> x)
                -> std::tuple<double, std::vector<double>> {
            // unpack density
            for (size_t i=0; i<x.size(); ++i) electron_density(i) = x[i];
            electron_density *= electron_density;
            double Nbar = integrate(electron_density, box);
            electron_density *= (N/Nbar);
            // get functional potentials
            double f;
            Double3D pot(electron_density);
            std::tie(f, pot) = energy_potential(false);
            // compute gradient
            Double3D grad = pot;
            pot *= electron_density;
            double mu = 1.0 / N * integrate(pot, box);
            grad -= mu;
            for (size_t i=0; i<x.size(); ++i) grad(i) *= x[i];
            grad *= (2.0 * N / Nbar * box.volume() / electron_density.size());
            // pack gradient and return
            std::vector<double> g(x);
            for (size_t i=0; i<x.size(); ++i) g[i] = grad(i);
            return {f, g};
        };

        bool converged;
        double f;
        size_t iterations;
        std::tie(converged, f, x, iterations) = 
            two_point_steepest_descent(compute_f_and_g, x);

    return (*this);
}

}
