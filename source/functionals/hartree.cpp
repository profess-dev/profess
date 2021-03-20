// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "hartree.hpp"

#include "fourier.hpp"

namespace profess
{

Hartree::Hartree(Box box)
    : Functional(box)
{
}

std::string Hartree::name()
{
    return "hartree";
}

void Hartree::set_box(std::array<std::array<double,3>,3> box)
{
    _box = Box(box);
}

double Hartree::energy(Double3D den)
{
    return std::get<0>(energy_potential(den));
}

std::tuple<double,Double3D> Hartree::energy_potential(Double3D den)
{
    const double fourpi = 4.0 * M_PI;
    Complex3D pot_ft = fourier_transform(den);
    Double3D k = wave_numbers(pot_ft.shape(), _box);
    // set k=0 component to zero. appropriate if the density integrates
    // to zero over the box (for example, if neutralized by a uniform 
    // background charge).
    pot_ft(0) = 0.0;
    pot_ft.set_elements(
        [fourpi, &pot_ft, &k](size_t i) {
            return fourpi * pot_ft(i) / (k(i)*k(i));
        },
        1,
        pot_ft.size());
    Double3D pot = inverse_fourier_transform(pot_ft, den.shape());
    double ene = 0.5 * (den*=pot).sum() * _box.volume() / den.size();
    return {ene, pot};
}

std::array<std::array<double,3>,3> Hartree::stress(Double3D den)
{
    std::array<std::array<double,3>,3> stress =
            {{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
    const double fourpi = 4.0 * M_PI;
    const double eightpi = 8.0 * M_PI;
    double hartree_energy = energy(den);
    Complex3D den_ft = fourier_transform(den);
    Double3D wvx = wave_vectors_x(den_ft.shape(), _box);
    Double3D wvy = wave_vectors_y(den_ft.shape(), _box);
    Double3D wvz = wave_vectors_z(den_ft.shape(), _box);
    Double3D wn = wave_numbers(den_ft.shape(), _box);
    for (size_t u=0; u<den_ft.shape()[0]; ++u) {
        for (size_t v=0; v<den_ft.shape()[1]; ++v) {
            for (size_t w=0; w<den_ft.shape()[2]; ++w) {
                // skip k=0 component. appropriate if density
                // integrates to zero over the box (for example,
                // if neutralized by a uniform background charge).
                if (u==0 and v==0 and w==0) continue;
                const double kx = wvx(u,v,w);
                const double ky = wvy(u,v,w);
                const double kz = wvz(u,v,w);
                const double kk = wn(u,v,w) * wn(u,v,w);
                const std::complex<double> ftden = den_ft(u,v,w);
                // most of the time, an extra factor of two is necessary
                // because only half of the complex fourier components are
                // stored along the third axis. the two exceptions are the
                // k3=0 plane (always) and the outermost k3 plane (if the
                // the real-space array size is even along that dimension).
                double factor;
                if (w==0 or (den.shape()[2]%2==0 and w+1==den_ft.shape()[2])) {
                    factor = fourpi;
                } else {
                    factor = eightpi;
                }
                const double factor_ftden2_over_k4 = 
                        factor * (ftden.real()*ftden.real() 
                                    + ftden.imag()*ftden.imag()) / (kk*kk);
                stress[0][0] += kx * kx * factor_ftden2_over_k4; // s_xx 
                stress[1][1] += ky * ky * factor_ftden2_over_k4; // s_yy
                stress[2][2] += kz * kz * factor_ftden2_over_k4; // s_zz
                stress[1][2] += ky * kz * factor_ftden2_over_k4; // s_yz
                stress[0][2] += kx * kz * factor_ftden2_over_k4; // s_xz
                stress[0][1] += kx * ky * factor_ftden2_over_k4; // s_xy
            }
        }
    }
    // finish stress tensor
    const double e_over_v = hartree_energy / _box.volume();
    for (size_t i=0; i<3; ++i) {
        stress[i][i] -= e_over_v;
        for (size_t j=i+1; j<3; ++j) {
            stress[j][i] = stress[i][j];
        }
    }
    return stress;
}

}
