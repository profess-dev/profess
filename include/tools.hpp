// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>

namespace profess
{

double interpolate_cubic_lagrange(std::vector<double> y, double x);
double interpolate_cubic_lagrange(
    std::vector<double> y, double x, double x0, double dx);

double interpolate_derivative_cubic_lagrange(std::vector<double> y, double x);
double interpolate_derivative_cubic_lagrange(
    std::vector<double> y, double x, double x0, double dx);

}
