// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "wang_govind_carter.hpp"

namespace profess
{

WangGovindCarter::WangGovindCarter(
        Box box,
        std::array<size_t,3> grid_shape,
        double den0,
        double alpha,
        double beta,
        double gamma)
    : Functional(box),
      _den0(den0),
      _alpha(alpha),
      _beta(beta),
      _gamma(gamma)
{
    _kedf_data = libkedf_init_();
    std::array<std::array<double,3>,3> vectors = box.vectors();
    int x = grid_shape[0];
    int y = grid_shape[1];
    int z = grid_shape[2];
    libkedf_initialize_grid_(
        _kedf_data, &x, &y, &z,
        vectors[0].data(), vectors[1].data(), vectors[2].data());
    libkedf_initialize_wgc2nd_custom_(
        _kedf_data, &_den0, &_alpha, &_beta, &_gamma);
}

std::string WangGovindCarter::name()
{
    return "WangGovindCarter";
}

void WangGovindCarter::set_box(std::array<std::array<double,3>,3> box)
{
    _box = Box(box);
}

double WangGovindCarter::energy(Double3D den)
{
    return std::get<0>(energy_potential(den));
}

std::tuple<double,Double3D> WangGovindCarter::energy_potential(Double3D den)
{
    double ene;
    Double3D pot(den.shape());
    libkedf_potential_(_kedf_data, den.data(), pot.data(), &ene);
    return {0.0, pot};
}

std::array<std::array<double,3>,3> WangGovindCarter::stress(Double3D den)
{
    std::array<std::array<double,3>,3> stress =
            {{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
    return stress;
}

}
