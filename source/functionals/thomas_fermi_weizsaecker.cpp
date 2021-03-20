// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "thomas_fermi_weizsaecker.hpp"

namespace profess
{

ThomasFermiWeizsaecker::ThomasFermiWeizsaecker(
        Box box,
        double a,
        double b,
        double tiny_den)
    : KineticGGA(
        box,
        [a,b,tiny_den](Double3D den, Double3D grad_dot_grad)
                -> std::tuple<Double3D,Double3D,Double3D> {
            const double c0 = 0.3*std::pow(3.0*M_PI*M_PI, 2.0/3.0);
            Double3D den23 = den;
            den23.compute_pow(2.0/3.0);
            Double3D integrand(den.shape());
            integrand.set_elements(
                [&](size_t i) {
                    return den(i)<tiny_den ? 0.0 :
                        a*c0*den(i)*den23(i)
                            + b*0.125*grad_dot_grad(i)/den(i);
                });
            Double3D deriv_1(den.shape());
            deriv_1.set_elements(
                [&](size_t i) {
                    return den(i)<tiny_den ? 0.0 :
                        a*5.0/3.0*c0*den23(i)
                            - b*0.125*grad_dot_grad(i)/(den(i)*den(i));
                });
            Double3D deriv_2(den.shape());
            deriv_2.set_elements(
                [&](size_t i) {
                    return den(i)<tiny_den ? 0.0 : b*0.125/den(i);
                });
            return {integrand, deriv_1, deriv_2};
        },
        [a,b,tiny_den](Double3D den, Double3D grad_dot_grad)
                -> std::tuple<Double3D,Double3D> {
            const double c0 = 0.3*std::pow(3.0*M_PI*M_PI, 2.0/3.0);
            Double3D den23 = den;
            den23.compute_pow(2.0/3.0);
            auto deriv_12 = Double3D(den.shape()).set_elements(
                [&](size_t i) {
                    return den(i)<tiny_den ? 0.0 :
                        - b*0.125/(den(i)*den(i));
                });
            auto deriv_22 = Double3D(den.shape()).fill(0.0);
            return {deriv_12, deriv_22};
        })
{
}

std::string ThomasFermiWeizsaecker::name()
{
    return "thomas-fermi-weizsaecker";
}

}
