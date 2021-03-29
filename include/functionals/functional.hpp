// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "array.hpp"
#include "box.hpp"

namespace profess
{

using Double3D = deft::Array<double>;
using Complex3D = deft::Array<std::complex<double>>;
using Box = deft::Box;

class Functional
{
public:
    Functional(Box box);
    virtual std::string name() = 0;
    virtual void set_box(std::array<std::array<double,3>,3> box) = 0;
    virtual double energy(Double3D den) = 0;
    virtual std::tuple<double, Double3D> energy_potential(Double3D den) = 0;
    virtual std::array<std::array<double,3>,3> stress(Double3D den) = 0;
protected:
    Box _box;
};

}
