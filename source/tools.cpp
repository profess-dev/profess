// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "tools.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace profess
{

double interpolate_cubic_lagrange(std::vector<double> y, double x)
{
// finds y(x) by interpolating with cubic lagrange polynomials.
//
// by default, assumes y = {y(i)} for i = {0, 1, ..., y.size()-1}.
//
// if, instead, y = {y(xi)} for xi = a + d * {0, 1, ..., y.size()-1},
// call the function with the signature:
//              y_interp = interpolate_cubic_lagrange(y, (x-a)/d).

    double y0, y1, y2, y3;             // yj values
    double x__x0, x__x1, x__x2, x__x3; // x-xj values

    if (x < 0.0) {
        throw std::runtime_error(
                "error: x<0 in 'interpolate_cubic_lagrange'.\n");
    }

    size_t i = static_cast<size_t>(std::floor(x));

    if (i == 0) {
        x__x0 = x;
        x__x1 = x - 1;
        x__x2 = x - 2;
        x__x3 = x - 3;
        y0 = y[0];
        y1 = y[1];
        y2 = y[2];
        y3 = y[3];
    } else if (i < y.size()-2) {
        double x__i = x - i; 
        x__x0 = 1 + x__i;
        x__x1 = x__i;
        x__x2 = x__i - 1;
        x__x3 = x__i - 2;
        y0 = y[i-1];
        y1 = y[i];
        y2 = y[i+1];
        y3 = y[i+2];
    } else if (i == y.size()-2) {
        double x__i = x - i; 
        x__x0 = 2 + x__i;
        x__x1 = 1 + x__i;
        x__x2 = x__i;
        x__x3 = x__i - 1;
        y0 = y[i-2];
        y1 = y[i-1];
        y2 = y[i];
        y3 = y[i+1];
    } else {
        throw std::runtime_error("error: in 'interpolate_cubic_lagrange', "
                                 "x outside of allowable range.\n");
    }

    constexpr double one_sixth = 1.0/6.0;
    constexpr double one_half = 0.5;
    return (- one_sixth * x__x1 * x__x2 * x__x3 * y0 
            + one_half  * x__x0 * x__x2 * x__x3 * y1
            - one_half  * x__x0 * x__x1 * x__x3 * y2
            + one_sixth * x__x0 * x__x1 * x__x2 * y3);
}

double interpolate_cubic_lagrange(
    std::vector<double> y, double x, double x0, double dx)
{
    return interpolate_cubic_lagrange(y, (x-x0)/dx);
}

double interpolate_derivative_cubic_lagrange(std::vector<double> y, double x)
{
// finds y'(x) by interpolating with cubic lagrange polynomials.
//
// by default, assumes y = {y(i)} for i = {0, 1, ..., y.size()-1}.
//
// if, instead, y = {y(xi)} for xi = a + d * {0, 1, ..., y.size()-1},
// call the function with the signature:
//          yp_interp = (1.0/d) * interpolate_cubic_lagrange(y, (x-a)/d).

    double y0, y1, y2, y3;             // yj values
    double x__x0, x__x1, x__x2, x__x3; // x-xj values

    if (x < 0.0) {
        throw std::runtime_error(
                "error: x<0 in 'interpolate_derivative_cubic_lagrange'.\n");
    }

    size_t i = static_cast<size_t>(std::floor(x));

    if (i == 0) {
        x__x0 = x;
        x__x1 = x - 1;
        x__x2 = x - 2;
        x__x3 = x - 3;
        y0 = y[0];
        y1 = y[1];
        y2 = y[2];
        y3 = y[3];
    } else if (i < y.size()-2) {
        double x__i = x - i; 
        x__x0 = 1 + x__i;
        x__x1 = x__i;
        x__x2 = x__i - 1;
        x__x3 = x__i - 2;
        y0 = y[i-1];
        y1 = y[i];
        y2 = y[i+1];
        y3 = y[i+2];
    } else if (i == y.size()-2) {
        double x__i = x - i; 
        x__x0 = 2 + x__i;
        x__x1 = 1 + x__i;
        x__x2 = x__i;
        x__x3 = x__i - 1;
        y0 = y[i-2];
        y1 = y[i-1];
        y2 = y[i];
        y3 = y[i+1];
    } else {
        throw std::runtime_error(
                "error: in 'interpolate_derivative_cubic_lagrange', "
                "x outside of allowable range.\n");
    }

    constexpr double one_sixth = 1.0/6.0;
    constexpr double one_half = 0.5;
    return (- one_sixth * (x__x2*x__x3 + x__x1*x__x3 + x__x1*x__x2) * y0
            + one_half  * (x__x2*x__x3 + x__x0*x__x3 + x__x0*x__x2) * y1
            - one_half  * (x__x1*x__x3 + x__x0*x__x3 + x__x0*x__x1) * y2
            + one_sixth * (x__x1*x__x2 + x__x0*x__x2 + x__x0*x__x1) * y3);
}

double interpolate_derivative_cubic_lagrange(
    std::vector<double> y, double x, double x0, double dx)
{
    return interpolate_derivative_cubic_lagrange(y, (x-x0)/dx) / dx;
}

}
