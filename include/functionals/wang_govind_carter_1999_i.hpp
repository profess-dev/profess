// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "kinetic_class_a.hpp"

namespace profess
{

class WangGovindCarter1999I : public KineticClassA
{
public:
    WangGovindCarter1999I(
        Box box,
        std::array<size_t,3> grid_shape,
        double den0,
        bool den0_immutable=false);
    std::string name();
private:
};

}
