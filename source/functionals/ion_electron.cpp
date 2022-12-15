// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "ion_electron.hpp"
#include "fourier.hpp"
#include "lattice_sum.hpp"
#include <algorithm>

namespace profess
{

IonElectron::IonElectron(
        std::array<size_t,3> grid_shape,
        Box box, Ions ions,
        int spline_order)
    : Functional::Functional(box),
      _potential(grid_shape),
      _ions(ions),
      _spline_order(spline_order)
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
    if (_spline_order>0) return _forces_spline(den);

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
    Complex3D str_fac(den_ft.shape());

    for (size_t i=0; i<_ions.count_types(); ++i) {

        const auto xyz_coords = _ions.xyz_coords_of_type(i);
        if (_spline_order<0) {
            str_fac = deft::structure_factor(den.shape(), _box, xyz_coords);
        } else {
            str_fac = deft::structure_factor_spline(
                den.shape(), _box, xyz_coords, _spline_order);
        }

        for (size_t u=0; u<den_ft.shape()[0]; ++u) {
        for (size_t v=0; v<den_ft.shape()[1]; ++v) {
        for (size_t w=0; w<den_ft.shape()[2]; ++w) {
    
            if (u==0 and v==0 and w==0) continue;

            const std::array<double,3> kv =
                {wvx(u,v,w), wvy(u,v,w), wvz(u,v,w)};
            const double k = wn(u,v,w);
            const double dv_dk = _ions.ft_potential_derivatives()[i](k);

            // third line of 't' expression requires explanation.
            // usually, factor of two is applied b/c of half-complex fft format.
            // but, must avoid double-counting for the edge cases noted.
            const double t = dv_dk / k
                    * (std::conj(den_ft(u,v,w)) * str_fac(u,v,w)).real()
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
            _potential.shape(),
            _box,
            loc,
            [this, i](double x, double y, double z) {
                return _ions.ft_potentials()[i](std::sqrt(x*x+y*y+z*z));
            },
            _spline_order);
    }
}

std::vector<std::array<double,3>> IonElectron::_forces_spline(Double3D den)
{
    std::vector<std::array<double,3>> forces(_ions.count(), {0.0,0.0,0.0});
    int N0 = den.shape()[0];
    int N1 = den.shape()[1];
    int N2 = den.shape()[2];
    Complex3D den_ft = fourier_transform(den);

    // TODO: move this eventually -----
    // get inverse of box vector matrix
    std::array<std::array<double,3>,3> ainv;
    std::array<std::array<double,3>,3> ainv_trans;
    auto b = _box.recip_vectors();
    for (size_t i=0; i<3; ++i) {
        for (size_t j=0; j<3; ++j) {
            // no transpose because b = [b1 b2 b3]^T
            ainv[i][j] = b[i][j] / (2.0*M_PI);
            ainv_trans[i][j] = b[j][i] / (2.0*M_PI);
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
    // compute fractional coordinates
    auto xyz_coords = _ions.xyz_coords();
    std::vector<std::array<double,3>> frac_coords(xyz_coords.size());
    for (size_t i=0; i<frac_coords.size(); ++i) {
        frac_coords[i] = m_dot_v(ainv, xyz_coords[i]);
    }
    // TODO: end delete eventually -----

    for (size_t i=0; i<_ions.count_types(); ++i) {

        Complex3D W_ft(den_ft.shape());
        auto wn = wave_numbers(W_ft.shape(), _box);
        for (int n0=0; n0<W_ft.shape()[0]; ++n0) {
            auto b0 = deft::exponential_spline_b(n0,N0,_spline_order);
            for (int n1=0; n1<W_ft.shape()[1]; ++n1) {
                auto b1 = deft::exponential_spline_b(n1,N1,_spline_order);
                for (int n2=0; n2<W_ft.shape()[2]; ++n2) {
                    auto b2 = deft::exponential_spline_b(n2,N2,_spline_order);
                    double k = wn(n0,n1,n2);
                    W_ft(n0,n1,n2) = den_ft(n0,n1,n2) / _box.volume();
                    W_ft(n0,n1,n2) *= N0*N1*static_cast<double>(N2)*b0*b1*b2;
                    W_ft(n0,n1,n2) *= _ions.ft_potentials()[i](k);
                    W_ft(n0,n1,n2) = std::conj(W_ft(n0,n1,n2));
                }
            }
        }
        auto W = inverse_fourier_transform(W_ft, den.shape());

        for (size_t a=0; a<frac_coords.size(); ++a) {

            if (_ions.type_ids()[a] != i) continue;

            double u0 = frac_coords[a][0]*N0;
            double u1 = frac_coords[a][1]*N1;
            double u2 = frac_coords[a][2]*N2;
            int floor0 = static_cast<int>(u0);
            int floor1 = static_cast<int>(u1);
            int floor2 = static_cast<int>(u2);
            auto M0 = deft::cardinal_b_spline_values(u0-floor0, _spline_order);
            auto M1 = deft::cardinal_b_spline_values(u1-floor1, _spline_order);
            auto M2 = deft::cardinal_b_spline_values(u2-floor2, _spline_order);
            auto M0p = deft::cardinal_b_spline_derivatives(u0-floor0, _spline_order);
            auto M1p = deft::cardinal_b_spline_derivatives(u1-floor1, _spline_order);
            auto M2p = deft::cardinal_b_spline_derivatives(u2-floor2, _spline_order);
            std::array<double,3> deriv_frac = {0.0,0.0,0.0};
            for (int i0=0; i0<_spline_order; ++i0) {
                int l0 = (i0-floor0)%N0 + (i0<floor0)*N0;
                for (int i1=0; i1<_spline_order; ++i1) {
                    int l1 = (i1-floor1)%N1 + (i1<floor1)*N1;
                    for (int i2=0; i2<_spline_order; ++i2) {
                        int l2 = (i2-floor2)%N2 + (i2<floor2)*N2;
                        deriv_frac[0] += N0*M0p[i0]*M1[i1]*M2[i2] * W(l0,l1,l2);
                        deriv_frac[1] += M0[i0]*N1*M1p[i1]*M2[i2] * W(l0,l1,l2);
                        deriv_frac[2] += M0[i0]*M1[i1]*N2*M2p[i2] * W(l0,l1,l2);
                    }
                }
            }
            forces[a] = m_dot_v(ainv_trans, deriv_frac);
            for (int j=0; j<3; j++)
                forces[a][j] = -forces[a][j] * _box.volume()/(N0*N1*N2);
        }
    }
    return forces;
}

}
