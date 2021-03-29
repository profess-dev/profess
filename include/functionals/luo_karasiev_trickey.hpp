// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "kinetic_gga.hpp"

namespace profess
{

class LuoKarasievTrickey : public KineticGGA
{
public:
    LuoKarasievTrickey(Box box, double a=1.3, double tiny_den=1e-12);
    std::string name();
};

}
