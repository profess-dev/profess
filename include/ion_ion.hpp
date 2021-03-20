// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "box.hpp"

#include "ions.hpp"

namespace profess
{

using Box = deft::Box;

class IonIon
{
public:
    std::string name();
    double energy(
            const Box box,
            const Ions& ions);
    std::vector<std::array<double,3>> forces(
            const Box box,
            const Ions& ions);
    std::array<std::array<double,3>,3> stress(
            const Box box,
            const Ions& ions);
};

}
