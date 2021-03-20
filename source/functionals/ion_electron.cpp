// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "ion_electron.hpp"
#include "fourier.hpp"
#include <algorithm>

namespace profess
{

IonElectron::IonElectron(std::array<size_t,3> grid_shape, Box box, Ions ions)
    : Functional::Functional(box),
      _potential(grid_shape),
      _ions(ions)
{
    _update_potential();
}

std::string IonElectron::name()
{
    return "ion_electron";
}

void IonElectron::set_box(std::array<std::array<double,3>,3> box)
{
    _box = Box(box);
    _update_potential();
}

double IonElectron::energy(Double3D den)
{
    return std::get<0>(energy_potential(den));
}

std::tuple<double, Double3D> IonElectron::energy_potential(
        Double3D density)
{
    density *= _potential;
    return {integrate(density,_box), _potential};
}

std::vector<std::array<double,3>> IonElectron::forces(Double3D den)
{
    std::vector<std::array<double,3>> forces(_ions.count(), {0.0,0.0,0.0});

    const size_t max_w = den.shape()[2]; 
    const bool max_w_even = (max_w%2 == 0);

    Complex3D den_ft = fourier_transform(den);
    Double3D wvx = wave_vectors_x(den_ft.shape(), _box);
    Double3D wvy = wave_vectors_y(den_ft.shape(), _box);
    Double3D wvz = wave_vectors_z(den_ft.shape(), _box);
    Double3D wn = wave_numbers(den_ft.shape(), _box);

    for (size_t i=0; i<_ions.count(); ++i) {

        const std::array<double,3> rv = _ions.xyz_coords()[i];

        for (size_t u=0; u<den_ft.shape()[0]; ++u) {
        for (size_t v=0; v<den_ft.shape()[1]; ++v) {
        for (size_t w=0; w<den_ft.shape()[2]; ++w) {
    
            if (u==0 and v==0 and w==0) continue;

            const std::array<double,3> kv = 
                {wvx(u,v,w), wvy(u,v,w), wvz(u,v,w)};
            const double k_dot_r = kv[0]*rv[0] + kv[1]*rv[1] + kv[2]*rv[2];
            const double k = wn(u,v,w);
            const double vk = _ions.ft_potentials()[_ions.type_ids()[i]](k);
            const std::complex<double> imag = {0.0,1.0};

            // third line of 't' expression requires explanation.
            // usually, factor of two applied b/c of half-complex fft format.
            // but, must avoid double-counting for some edge cases.
            const double t = vk
                    * (std::exp(imag*k_dot_r) * den_ft(u,v,w)).imag()
                    * ((w==0 or (max_w_even and w==max_w)) ? 1.0 : 2.0);

            for (size_t c=0; c<3; ++c) forces[i][c] += t * kv[c];
        }}}
    }

    return forces;
}

std::array<std::array<double,3>,3> IonElectron::stress(Double3D den)
{
    std::array<std::array<double,3>,3> stress =
            {{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};

    const size_t max_w = den.shape()[2]; 
    const bool max_w_even = (max_w%2 == 0);

    Complex3D den_ft = fourier_transform(den);
    Double3D wvx = wave_vectors_x(den_ft.shape(), _box);
    Double3D wvy = wave_vectors_y(den_ft.shape(), _box);
    Double3D wvz = wave_vectors_z(den_ft.shape(), _box);
    Double3D wn = wave_numbers(den_ft.shape(), _box);

    for (size_t i=0; i<_ions.count_types(); ++i) {

        const std::vector<std::array<double,3>> xyz_coords =
                _ions.xyz_coords_of_type(i);

        for (size_t u=0; u<den_ft.shape()[0]; ++u) {
        for (size_t v=0; v<den_ft.shape()[1]; ++v) {
        for (size_t w=0; w<den_ft.shape()[2]; ++w) {
    
            if (u==0 and v==0 and w==0) continue;

            const std::array<double,3> kv =
                {wvx(u,v,w), wvy(u,v,w), wvz(u,v,w)};
            const double k = wn(u,v,w);
            const double dv_dk = _ions.ft_potential_derivatives()[i](k);

            std::complex<double> structure_factor = 0.0;
            for (size_t j=0; j<xyz_coords.size(); ++j) {
                const std::complex<double> imag = {0.0,1.0};
                const std::array<double,3> xyz = xyz_coords[j];
                const double k_dot_r =
                        kv[0]*xyz[0] + kv[1]*xyz[1] + kv[2]*xyz[2];
                structure_factor += std::exp(-imag * k_dot_r);
            }

            // third line of 't' expression requires explanation.
            // usually, factor of two is applied b/c of half-complex fft format.
            // but, must avoid double-counting for the edge cases noted.
            const double t = dv_dk / k
                    * (std::conj(den_ft(u,v,w)) * structure_factor).real()
                    * ((w==0 or (max_w_even and w==max_w)) ? 1.0 : 2.0);
            for (size_t a=0; a<3; ++a) {
                for (size_t b=0; b<=a; ++b) {
                    stress[a][b] -= t * kv[a] * kv[b];
                }
            }
        }}}
    }

    // fill remaining off-diagonal elements by symmetry
    for (size_t a=0; a<3; ++a) {
        for (size_t b=a+1; b<3; ++b) {
            stress[a][b] = stress[b][a];
        }
    }

    // finalize stress tensor
    const double s = energy(den);
    for (size_t a=0; a<3; ++a) {
        stress[a][a] -= s;
        for (size_t b=0; b<3; ++b) {
            stress[a][b] /= _box.volume();
        }
    }
    return stress;
}

void IonElectron::_update_potential()
{
    _potential.fill(0.0);
    for (size_t i=0; i<_ions.count_types(); ++i) {

        // extract locations of ions of type i
        size_t num_i = _ions.count_of_type(i);
        std::vector<std::array<double,3>> loc(num_i);
        size_t j = 0;
        for (size_t ii=0; ii<_ions.count(); ++ii) {
            if (_ions.type_ids()[ii] == i) {
                loc[j] = _ions.xyz_coords()[ii];
                j += 1;
            }
        }
        // get potential from ions of type i
        _potential += array_from_lattice_sum(
            _potential.shape(), _box, loc, _ions.ft_potentials()[i]);
    }
}

}
