// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include <cmath>

#include "c_real_space_electrostatic_sum.h"
#include "ion_ion.hpp"

namespace profess
{

std::string IonIon::name()
{
    return "ion_ion";
}

double IonIon::energy(
        const Box box,
        const Ions& ions)
{
    // compute max distance between lattice planes
    const double d0 = 2.0*M_PI/box.recip_lengths()[0];
    const double d1 = 2.0*M_PI/box.recip_lengths()[1];
    const double d2 = 2.0*M_PI/box.recip_lengths()[2];
    const double h_max = std::max(d0, std::max(d1, d2));

    const double r_c = 250.0;
    const double r_d = std::sqrt(1.0/3.0*h_max*r_c);

    const int num_ions = ions.count();
    std::vector<double> rx(num_ions), ry(num_ions), rz(num_ions);
    std::vector<double> ion_charges(num_ions);
    for (int i=0; i<num_ions; ++i) {
        const std::array<double,3> r = ions.xyz_coords()[i];
        rx[i] = r[0];
        ry[i] = r[1];
        rz[i] = r[2];
        ion_charges[i] = ions.charges()[ions.type_ids()[i]];
    }

    double ene;
    c_real_space_electrostatic_sum_energy(
            box.vectors()[0].data(),
            box.vectors()[1].data(),
            box.vectors()[2].data(),
            &num_ions,
            rx.data(), ry.data(), rz.data(),
            ion_charges.data(),
            &r_c,
            &r_d,
            &ene);

    return ene;
}

std::vector<std::array<double,3>> IonIon::forces(
        const Box box,
        const Ions& ions)
{
    // compute max distance between lattice planes
    const double d0 = 2.0*M_PI/box.recip_lengths()[0];
    const double d1 = 2.0*M_PI/box.recip_lengths()[1];
    const double d2 = 2.0*M_PI/box.recip_lengths()[2];
    const double h_max = std::max(d0, std::max(d1, d2));

    const double r_c = 250.0;
    const double r_d = std::sqrt(1.0/3.0*h_max*r_c);

    const int num_ions = ions.count();
    std::vector<double> rx(num_ions), ry(num_ions), rz(num_ions);
    std::vector<double> ion_charges(num_ions);
    for (int i=0; i<num_ions; ++i) {
        const std::array<double,3> r = ions.xyz_coords()[i];
        rx[i] = r[0];
        ry[i] = r[1];
        rz[i] = r[2];
        ion_charges[i] = ions.charges()[ions.type_ids()[i]];
    }

    std::vector<double> fx(num_ions), fy(num_ions), fz(num_ions);

    c_real_space_electrostatic_sum_force(
            box.vectors()[0].data(),
            box.vectors()[1].data(),
            box.vectors()[2].data(),
            &num_ions,
            rx.data(), ry.data(), rz.data(),
            ion_charges.data(),
            &r_c,
            &r_d,
            fx.data(), fy.data(), fz.data());

    std::vector<std::array<double,3>> forces(num_ions);
    for (int i=0; i<num_ions; ++i)
        forces[i] = {fx[i], fy[i], fz[i]};

    return forces;
}

std::array<std::array<double,3>,3> IonIon::stress(
        const Box box,
        const Ions& ions)
{
    // compute max distance between lattice planes
    const double d0 = 2.0*M_PI/box.recip_lengths()[0];
    const double d1 = 2.0*M_PI/box.recip_lengths()[1];
    const double d2 = 2.0*M_PI/box.recip_lengths()[2];
    const double h_max = std::max(d0, std::max(d1, d2));

    const double r_c = 250.0;
    const double r_d = std::sqrt(1.0/3.0*h_max*r_c);

    std::array<std::array<double,3>,3> stress =
            {{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};

    const int num_ions = ions.count();
    std::vector<double> rx(num_ions), ry(num_ions), rz(num_ions);
    std::vector<double> ion_charges(num_ions);
    for (int i=0; i<num_ions; ++i) {
        const std::array<double,3> r = ions.xyz_coords()[i];
        rx[i] = r[0];
        ry[i] = r[1];
        rz[i] = r[2];
        ion_charges[i] = ions.charges()[ions.type_ids()[i]];
    }

    std::array<double,6> stress_voight;

    c_real_space_electrostatic_sum_stress(
            box.vectors()[0].data(),
            box.vectors()[1].data(),
            box.vectors()[2].data(),
            &num_ions,
            rx.data(), ry.data(), rz.data(),
            ion_charges.data(),
            &r_c,
            &r_d,
            stress_voight.data());

    stress[0][0] = stress_voight[0];
    stress[1][1] = stress_voight[1];
    stress[2][2] = stress_voight[2];
    stress[1][2] = stress[2][1] = stress_voight[3];
    stress[0][2] = stress[2][0] = stress_voight[4];
    stress[0][1] = stress[1][0] = stress_voight[5];

    return stress;
}

}
