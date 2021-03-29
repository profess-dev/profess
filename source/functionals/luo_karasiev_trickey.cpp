// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "luo_karasiev_trickey.hpp"

namespace profess
{

LuoKarasievTrickey::LuoKarasievTrickey(Box box, double a, double tiny_den)
    : KineticGGA(
        box,
        [a,tiny_den](Double3D den, Double3D grad_dot_grad)
                -> std::tuple<Double3D,Double3D,Double3D> {
            // compute useful quantities
            const double c0 = 0.3*std::pow(3.0*M_PI*M_PI, 2.0/3.0);
            const double cs = 2.0*std::pow(3.0*M_PI*M_PI, 1.0/3.0);
            Double3D den23 = den;
            den23.compute_pow(2.0/3.0);
            Double3D s = (grad_dot_grad/(cs*cs*den*den*den23)).compute_sqrt();
            auto cosh = Double3D(den.shape()).set_elements(
                [a,&s](size_t i) { return std::cosh(a*s(i)); });
            auto tanh = Double3D(den.shape()).set_elements(
                [a,&s](size_t i) { return std::tanh(a*s(i)); });
            // compute integrand
            auto integrand = Double3D(den.shape()).set_elements(
                [&](size_t i) {
                    return (den(i) < tiny_den) ? 0.0 :
                        grad_dot_grad(i)/(8.0*den(i))
                            + c0*den(i)*den23(i)/cosh(i);
                });
            // compute derivative wrt den
            auto deriv_1 = Double3D(den.shape()).set_elements(
                [&](size_t i) {
                    return (den(i) < tiny_den) ? 0.0 :
                        -grad_dot_grad(i)/(8.0*den(i)*den(i))
                            + c0*den23(i)/cosh(i)*(
                                5.0/3.0+4.0/3.0*a*s(i)*tanh(i));
                });
            // compute derivative wrt grad_dot_grad
            auto deriv_2 = Double3D(den.shape()).set_elements(
                [&](size_t i) {
                    if (den(i) < tiny_den) return 0.0;
                    const double tanh_over_s = (s(i) < 1e-6) ? 
                        a - a*a*a*s(i)*s(i)/3.0 : // series for small s
                        std::tanh(a*s(i))/s(i);   // otherwise tanh(as)/s
                    return 1.0/(8.0*den(i))
                        - 3.0/(80.0*den(i))*a/cosh(i)*tanh_over_s;

                });
            return {integrand, deriv_1, deriv_2};
        },
        [a,tiny_den](Double3D den, Double3D grad_dot_grad)
                -> std::tuple<Double3D,Double3D> {
            // compute useful quantities
            const double cs = 2.0*std::pow(3.0*M_PI*M_PI, 1.0/3.0);
            Double3D den23 = den;
            den23.compute_pow(2.0/3.0);
            Double3D s = (grad_dot_grad/(cs*cs*den*den*den23)).compute_sqrt();
            auto cosh = Double3D(den.shape()).set_elements(
                [a,&s](size_t i) { return std::cosh(a*s(i)); });
            auto tanh = Double3D(den.shape()).set_elements(
                [a,&s](size_t i) { return std::tanh(a*s(i)); });
            // compute deriv_12
            auto deriv_12 = Double3D(den.shape()).set_elements(
                [&](size_t i) {
                    if (den(i) < tiny_den) return 0.0;
                    const double tanh_over_s = (s(i) < 1e-6) ? 
                        a - a*a*a*s(i)*s(i)/3.0 : // series for small s
                        std::tanh(a*s(i))/s(i);   // otherwise tanh(as)/s
                    return - 0.125/(den(i)*den(i))
                        - a/(80.0*den(i)*den(i)*cosh(i))*tanh_over_s
                        + 4.0*a*a/(80.0*den(i)*den(i)*cosh(i))*(
                            2.0/(cosh(i)*cosh(i)) - 1.0);
                });
            auto deriv_22 = Double3D(den.shape()).set_elements(
                [&](size_t i) {
                    if (den(i) < tiny_den) return 0.0;
                    const double d = (s(i) < 1e-6) ? 
                        -5.0/3.0+6.0/5.0*a*a*s(i)*s(i) : // small s
                        (2.0/(cosh(i)*cosh(i))-1.0-tanh(i)/(a*s(i))) /
                            (a*a*s(i)*s(i)); // otherwise
                    return -3.0*a*a*a*a / (
                        160.0*cs*cs*den(i)*den(i)*den(i)*den23(i)*cosh(i))
                        * d;

                });
            return {deriv_12, deriv_22};
        })
{
}

std::string LuoKarasievTrickey::name()
{
    return "luo-karasiev-trickey";
}

}
