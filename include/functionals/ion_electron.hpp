// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "functional.hpp"
#include "ions.hpp"

namespace profess
{

class IonElectron : public Functional
{
public:
    IonElectron(
        std::array<size_t,3> grid_shape,
        Box box,
        Ions ions,
        int spline_order=-1);
    std::string name() override;
    void set_box(std::array<std::array<double,3>,3> box) override;
    double energy(Double3D den) override;
    std::tuple<double, Double3D> energy_potential(Double3D den) override;
    std::vector<std::array<double,3>> forces(Double3D den);
    std::array<std::array<double,3>,3> stress(Double3D den) override;
private:
    Double3D _potential;
    Ions _ions;
    void _update_potential();
    int _spline_order;
    std::vector<std::array<double,3>> _forces_spline(Double3D den);
};

}
