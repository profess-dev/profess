// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "kinetic_gga.hpp"

namespace profess
{

class ThomasFermiWeizsaecker : public KineticGGA
{
public:
    ThomasFermiWeizsaecker(
        Box box,
        double a=1.0,
        double b=1.0,
        double tiny_den=1e-12);
    std::string name();
};

}
