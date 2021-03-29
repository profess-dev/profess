// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "weizsaecker.hpp"

#include "fourier.hpp"

namespace profess
{

std::string Weizsaecker::name()
{
    return "weizsaecker";
}

void Weizsaecker::set_box(std::array<std::array<double,3>,3> box)
{
    _box = Box(box);
}

double Weizsaecker::energy(Double3D den)
{
    return std::get<0>(energy_potential(den));
}

std::tuple<double, Double3D> Weizsaecker::energy_potential(Double3D den)
{
    Double3D pot(den);
    pot.compute_sqrt();
    pot = -0.5 * pot * laplacian(pot, _box);
    const double ene = pot.sum() * _box.volume() / pot.size();
    pot /= den;
    return {ene, pot};
}

std::array<std::array<double,3>,3> Weizsaecker::stress(Double3D den)
{
    std::array<std::array<double,3>,3> stress =
            {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    // store fourier transform of sqrt(den)
    den.compute_sqrt();
    Complex3D sqrt_ft = fourier_transform(den);
    // get grad_x, grad_y, and grad_z
    Double3D grad_x =
        inverse_fourier_transform(gradient_x(sqrt_ft,_box), den.shape());
    Double3D grad_y =
        inverse_fourier_transform(gradient_y(sqrt_ft,_box), den.shape());
    Double3D grad_z =
        inverse_fourier_transform(gradient_z(sqrt_ft,_box), den.shape());
    // s_xx
    den = grad_x * grad_x;
    stress[0][0] = -1.0/_box.volume()*(den.sum()*_box.volume()/den.size());
    // s_yy
    den = grad_y * grad_y;
    stress[1][1] = -1.0/_box.volume()*(den.sum()*_box.volume()/den.size());
    // s_zz
    den = grad_z * grad_z;
    stress[2][2] = -1.0/_box.volume()*(den.sum()*_box.volume()/den.size());
    // s_yz
    den = grad_y * grad_z;
    stress[1][2] = -1.0/_box.volume()*(den.sum()*_box.volume()/den.size());
    stress[2][1] = stress[1][2];
    // s_xz
    den = grad_x * grad_z;
    stress[0][2] = -1.0/_box.volume()*(den.sum()*_box.volume()/den.size());
    stress[2][0] = stress[0][2];
    // s_xy
    den = grad_x * grad_y;
    stress[0][1] = -1.0/_box.volume()*(den.sum()*_box.volume()/den.size());
    stress[1][0] = stress[0][1];
    return stress;
}

}
