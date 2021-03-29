// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "functional.hpp"

namespace profess
{

class Weizsaecker : public Functional
{
public:
    using Functional::Functional; // inherit base class constructor
    std::string name() override;
    void set_box(std::array<std::array<double,3>,3> cell) override;
    double energy(Double3D den) override;
    std::tuple<double, Double3D> energy_potential(Double3D den) override;
    std::array<std::array<double,3>,3> stress(Double3D den) override;
};

}
