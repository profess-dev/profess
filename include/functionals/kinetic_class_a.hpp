// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "functional.hpp"

namespace profess
{

class KineticClassA : public Functional
{
public:
    KineticClassA(
        Box box,
        std::array<size_t,3> grid_shape,
        double a,
        double b,
        std::function<double(double)> f,
        std::function<double(double)> fp,
        double den0,
        bool den0_immutable=false);
    std::string name() override;
    void set_box(std::array<std::array<double,3>,3> box) override;
    double energy(Double3D den) override;
    std::tuple<double, Double3D> energy_potential(Double3D den) override;
    std::array<std::array<double,3>,3> stress(Double3D den) override;
private:
    std::array<size_t,3> _grid_shape;
    const double _a;
    const double _b;
    const std::function<double(double)> _f;
    const std::function<double(double)> _fp;
    double _den0;
    const bool _den0_immutable;
    Double3D _kernel;
    void _generate_kernel();
    double _lindhard(double x);
    double _lindhard_deriv(double x);
};

}
