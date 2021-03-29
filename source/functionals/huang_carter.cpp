// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "huang_carter.hpp"

namespace profess
{

HuangCarter::HuangCarter(
        Box box,
        std::array<size_t,3> grid_shape,
        double den0,
        double lambda,
        double alpha,
        double beta,
        double refRatio)
    : Functional(box),
      _den0(den0),
      _lambda(lambda),
      _alpha(alpha),
      _beta(beta),
      _refRatio(refRatio)
{
    _kedf_data = libkedf_init_();
    std::array<std::array<double,3>,3> vectors = box.vectors();
    int x = grid_shape[0];
    int y = grid_shape[1];
    int z = grid_shape[2];
    libkedf_initialize_grid_(
        _kedf_data, &x, &y, &z,
        vectors[0].data(), vectors[1].data(), vectors[2].data());
    libkedf_initialize_hc_custom_(
        _kedf_data, &_den0, &_lambda, &_alpha, &_beta, &_refRatio);
}

std::string HuangCarter::name()
{
    return "HuangCarter";
}

void HuangCarter::set_box(std::array<std::array<double,3>,3> box)
{
    _box = Box(box);
}

double HuangCarter::energy(Double3D den)
{
    return std::get<0>(energy_potential(den));
}

std::tuple<double,Double3D> HuangCarter::energy_potential(Double3D den)
{
    double ene;
    Double3D pot(den.shape());
    libkedf_potential_(_kedf_data, den.data(), pot.data(), &ene);
    return {0.0, pot};
}

std::array<std::array<double,3>,3> HuangCarter::stress(Double3D den)
{
    std::array<std::array<double,3>,3> stress =
            {{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
    return stress;
}

}
