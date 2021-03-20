// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include <algorithm>

#include "minimizers.hpp"

namespace profess
{

std::tuple<
    bool,                // converged?
    double,              // f
    std::vector<double>, // x
    size_t>              // iterations
two_point_steepest_descent(
    std::function<std::tuple<double,std::vector<double>>(
                      std::vector<double>)> compute_f_and_g,
    std::vector<double> x,
    double f_tol,
    size_t f_window,
    double g_tol,
    size_t max_iter)
{
    // for initial x (x0), get function value and gradient (f0/g0)
    double f;
    std::vector<double> g(x.size());
    std::tie(f, g) = compute_f_and_g(x);
    // initialize vector for storing f values
    std::vector<double> f_history(f_window);
    f_history[0] = f;
    for (size_t j=1; j<f_window; ++j) {
        f_history[j] = f_history[j-1] + 10*f_tol;
    }
    // get x1/f1/g1
    double step = 1.0;
    auto x_prev = x;
    auto g_prev = g;
    for (size_t j=0; j<x.size(); ++j) x[j] = x[j] - step*g[j];
    std::tie(f, g) = compute_f_and_g(x);
    // loop until max_iter or convergence
    size_t iterations = 1;
    bool converged = false;
    for (size_t i=2; i<=max_iter; ++i) {
        iterations++;
        // compute step
        double dx_dot_dg = 0.0;
        double dg_dot_dg = 0.0;
        for (size_t j=0; j<x.size(); ++j) {
            double dgj = g[j] - g_prev[j];
            dx_dot_dg += (x[j]-x_prev[j])*dgj;
            dg_dot_dg += dgj*dgj;
        }
        const double step = dx_dot_dg / dg_dot_dg;
        // get xi/fi/gi
        x_prev = x;
        g_prev = g;
        for (size_t j=0; j<x.size(); ++j) x[j] = x[j] - step*g[j];
        std::tie(f,g) = compute_f_and_g(x);
        // check for convergence
        bool f_converged = std::all_of(
                f_history.begin(),
                f_history.end(),
                [f,f_tol](double e) { return std::abs(e-f) < f_tol; });
        bool g_converged = false;
        if (f_converged) {
            double g_dot_g = 0.0;
            for (size_t j=0; j<x.size(); ++j) g_dot_g += g[j]*g[j];
            if (g_dot_g < g_tol*g_tol) g_converged = true;
        }
        if (f_converged and g_converged) {
            converged = true;
            break;
        } else {
            for (size_t j=f_window-1; j>0; --j) {
                f_history[j] = f_history[j-1];
            }
            f_history[0] = f;
        }
    }
    return {converged, f, x, iterations};
}

}
