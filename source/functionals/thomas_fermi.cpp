// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "thomas_fermi.hpp"

namespace profess
{

std::string ThomasFermi::name()
{
    return "thomas_fermi";
}

void ThomasFermi::set_box(std::array<std::array<double,3>,3> box)
{
    _box = Box(box);
}

double ThomasFermi::energy(Double3D den)
{
    return std::get<0>(energy_potential(den));
}

std::tuple<double,Double3D> ThomasFermi::energy_potential(Double3D den)
{
    const double c0 = 0.3 * std::pow(3.0*M_PI*M_PI, 2.0/3.0);
    const double c0_53 = c0 * 5.0/3.0;
    Double3D pot(den);
    pot.set_elements(
        [c0_53,&pot](size_t i) {
            const double cbrt = std::cbrt(pot(i));
            return c0_53 * cbrt * cbrt;
        });
    return {(3.0/5.0)*(den*=pot).sum()*_box.volume()/den.size(), pot};
}

std::array<std::array<double,3>,3> ThomasFermi::stress(Double3D den)
{
    std::array<std::array<double,3>,3> stress =
            {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    double s = -2.0/(3.0*_box.volume()) * energy(den);
    for (int i=0; i<3; ++i) stress[i][i] = s;
    return stress;
}

}
