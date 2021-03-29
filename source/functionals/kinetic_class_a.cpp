// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "kinetic_class_a.hpp"
#include "thomas_fermi.hpp"
#include "weizsaecker.hpp"

#include "fourier.hpp"

namespace profess
{

KineticClassA::KineticClassA(
        Box box,
        std::array<size_t,3> grid_shape,
        double a,
        double b,
        std::function<double(double)> f,
        std::function<double(double)> fp,
        double den0,
        bool den0_immutable)
    : Functional::Functional(box),
      _grid_shape(grid_shape),
      _a(a),
      _b(b),
      _f(f),
      _fp(fp),
      _den0(den0),
      _den0_immutable(den0_immutable),
      _kernel({grid_shape[0], grid_shape[1], grid_shape[2]/2+1})
{
    _generate_kernel();
}

std::string KineticClassA::name()
{
    return "KineticClassA";
}

void KineticClassA::set_box(std::array<std::array<double,3>,3> box)
{
    if (_den0_immutable) {
        _box = Box(box);
    } else {
        _den0 *= _box.volume();  // multiply by old volume
        _box = Box(box);       // set new box
        _den0 /= _box.volume();  // divide by new volume
    }
    _generate_kernel();
}

double KineticClassA::energy(Double3D den)
{
    return std::get<0>(energy_potential(den));
}

std::tuple<double, Double3D> KineticClassA::energy_potential(Double3D den)
{
    // weizsaecker energy and potential
    Weizsaecker f_w(_box);
    double ene_w;
    Double3D pot_w(den);
    std::tie(ene_w, pot_w) = f_w.energy_potential(den);

    // thomas-fermi energy and potential
    ThomasFermi f_tf(_box);
    double ene_tf;
    Double3D pot_tf(den);
    std::tie(ene_tf, pot_tf) = f_tf.energy_potential(den);

    // compute X[n]
    Double3D pot(den);
    pot.compute_pow(_b);
    Complex3D pot_ft = fourier_transform(pot);
    pot_ft.set_elements(
        [&pot_ft, this](size_t i) { return pot_ft(i) * _kernel(i); });
    pot = inverse_fourier_transform(pot_ft, pot.shape());
    Double3D tmp(den);
    tmp.compute_pow(_a);
    tmp *= pot;
    const double X = integrate(tmp,_box) / ene_tf;

    tmp = den;
    tmp.compute_pow(_a-1.0);
    pot *= _a*tmp;

    tmp = den;
    tmp.compute_pow(_a);
    Complex3D tmp_ft = fourier_transform(tmp);
    tmp_ft.set_elements(
        [&tmp_ft, this](size_t i) { return tmp_ft(i) * _kernel(i); });
    tmp = inverse_fourier_transform(tmp_ft, tmp.shape());
    den.compute_pow(_b-1.0);
    pot += _b*den*tmp;

    // finalize energy and potential
    double ene = ene_w + ene_tf*_f(X);
    pot = pot_w + pot_tf*(_f(X)-_fp(X)*X) + _fp(X)*pot;
    return {ene, pot};
}

std::array<std::array<double,3>,3> KineticClassA::stress(Double3D den)
{
    // weizsaecker stress
    std::array<std::array<double,3>,3> stress = Weizsaecker(_box).stress(den);

    // thomas-fermi energy and stress, used later
    ThomasFermi f_tf(_box);
    double ene_tf = f_tf.energy(den);
    std::array<std::array<double,3>,3> stress_tf = f_tf.stress(den);

    // compute free-electron-gas parameters
    const double k0 = std::cbrt(3.0*M_PI*M_PI*_den0);

    // compute X[n]
    const double kernel_prefactor = M_PI*M_PI/k0
            * std::pow(_den0, 2.0-_a-_b) / (2.0*_a*_b*_fp(0.0));
    Double3D wrk(den);
    wrk.compute_pow(_b);
    Complex3D wrk_ft = fourier_transform(wrk);
    Double3D wave_vectors_k = wave_numbers(wrk_ft.shape(), _box);
    wrk_ft.set_elements(
        [&wave_vectors_k, k0, &wrk_ft, kernel_prefactor, this](size_t i) {
            const double eta = wave_vectors_k(i) / (2.0*k0);
            return wrk_ft(i) * kernel_prefactor
                * (1.0/_lindhard(eta)-3.0*eta*eta-1.0);
            });
    wrk = inverse_fourier_transform(wrk_ft, wrk.shape());
    Double3D tmp(den);
    tmp.compute_pow(_a);
    tmp *= wrk;
    const double X = integrate(tmp, _box) / ene_tf;

    // more contributions to stress tensor
    const double c = _den0_immutable ? (1.0-_a-_b) : (-2.0/3.0);
    for (size_t i=0; i<3; ++i) {
        stress[i][i] += (stress_tf[i][i]*_f(X)
                           - stress_tf[i][i]*_fp(X)*X
                           + c/_box.volume()*ene_tf*_fp(X)*X);
    }

    // loop to generate final contributions to stress tensor
    double s_xx=0.0, s_yy=0.0, s_zz=0.0, s_yz=0.0, s_xz=0.0, s_xy=0.0;
    wrk = den;
    wrk.compute_pow(_a);
    wrk_ft = fourier_transform(wrk);
    tmp = den;
    tmp.compute_pow(_b);
    Complex3D tmp_ft = fourier_transform(tmp);
    Double3D wave_vectors_kx = wave_vectors_x(wrk_ft.shape(), _box);
    Double3D wave_vectors_ky = wave_vectors_y(wrk_ft.shape(), _box);
    Double3D wave_vectors_kz = wave_vectors_z(wrk_ft.shape(), _box);
    for (size_t u=0; u<wrk_ft.shape()[0]; ++u) {
        for (size_t v=0; v<wrk_ft.shape()[1]; ++v) {
            for (size_t w=0; w<wrk_ft.shape()[2]; ++w) {
                // skip k=0 to avoid nan's (allowed bc the kernel vanishes)
                if (u==0 and v==0 and w==0) continue;
                const double k = wave_vectors_k(u,v,w);
                const double eta = k / (2.0*k0);
                const double lind = _lindhard(eta);
                const double d_lind = _lindhard_deriv(eta);
                // get   \tilde{n^a}_{-k} * \tilde{n^b}_{k}
                //     + \tilde{n^a}_{k}  * \tilde{n^b}_{-k}
                const double dena_denb =
                        2.0 * (tmp_ft(u,v,w).real()*wrk_ft(u,v,w).real()
                                + tmp_ft(u,v,w).imag()*wrk_ft(u,v,w).imag());
                double value = (-d_lind/(lind*lind)-6.0*eta)*dena_denb;
                // most of the time, the full dena_denb is
                // kept because of the half-complex storage format
                // along the third axis. the two exceptions are the
                // k3=0 plane and the outermost k3 plane (if the
                // the real-space array size is even along that
                // dimension), where double counting must be removed.
                if (w==0 or (wrk.shape()[2]%2==0 and w+1==wrk_ft.shape()[2]))
                    value *= 0.5;
                const double kx = wave_vectors_kx(u,v,w);
                const double ky = wave_vectors_ky(u,v,w);
                const double kz = wave_vectors_kz(u,v,w);
                const double extra = _den0_immutable ? 0.0 : eta/3.0;
                s_xx += value * (-kx*kx/(2.0*k*k0) + extra); 
                s_yy += value * (-ky*ky/(2.0*k*k0) + extra);
                s_zz += value * (-kz*kz/(2.0*k*k0) + extra);
                s_yz += value * (-ky*kz/(2.0*k*k0));
                s_xz += value * (-kx*kz/(2.0*k*k0));
                s_xy += value * (-kx*ky/(2.0*k*k0));
            }
        }
    }
    stress[0][0] += _fp(X) * kernel_prefactor * s_xx;
    stress[1][1] += _fp(X) * kernel_prefactor * s_yy;
    stress[2][2] += _fp(X) * kernel_prefactor * s_zz;
    stress[1][2] += _fp(X) * kernel_prefactor * s_yz;
    stress[0][2] += _fp(X) * kernel_prefactor * s_xz;
    stress[0][1] += _fp(X) * kernel_prefactor * s_xy;
    stress[2][1] = stress[1][2];
    stress[2][0] = stress[0][2];
    stress[1][0] = stress[0][1];

    return stress;
}

void KineticClassA::_generate_kernel()
{
    Double3D wave_vectors_k = wave_numbers(_kernel.shape(), _box);
    const double k0 = std::cbrt(3.0*M_PI*M_PI*_den0);
    const double factor = M_PI*M_PI/k0
            * std::pow(_den0,2.0-_a-_b) / (2.0*_a*_b*_fp(0.0));
    _kernel.set_elements(
        [&wave_vectors_k, k0, factor, this](size_t i) {
            const double eta = wave_vectors_k(i) / (2.0*k0);
            return factor * (1.0/_lindhard(eta) - 3.0*eta*eta - 1.0);
        });
}

double KineticClassA::_lindhard(double x)
{
    if (x < 1e-8) {
        return 1.0;
    } else {
        return 0.5 + (1.0-x*x)/(4.0*x)*std::log(std::abs((1.0+x)/(1.0-x)));
    }
}

double KineticClassA::_lindhard_deriv(double x)
{
    if (x < 1e-8) {
        return 0.0;
    } else {
        return 0.5/x-(1.0+x*x)/(4.0*x*x)*std::log(std::abs((1.0+x)/(1.0-x)));
    }
}

}
