// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include <functional>
#include <vector>

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
    double f_tol=1e-7,
    size_t f_window=2,
    double g_tol=1e-5,
    size_t max_iter=100);

}
