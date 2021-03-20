// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "perdew_zunger.hpp"

namespace profess
{

PerdewZunger::PerdewZunger(Box box)
    : LibXC(box, {1,10})
{
}

std::string PerdewZunger::name()
{
    return "pbe";
}

}
