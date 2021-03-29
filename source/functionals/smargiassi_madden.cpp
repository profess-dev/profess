// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "smargiassi_madden.hpp"

namespace profess
{

SmargiassiMadden::SmargiassiMadden(
        Box box,
        std::array<size_t,3> grid_shape,
        double den0,
        bool den0_immutable)
    : KineticClassA(
        box,
        grid_shape,
        0.5,
        0.5,
        [](double x){return 1.0+x;},
        [](double x){return 1.0;},
        den0,
        den0_immutable)
{
}

std::string SmargiassiMadden::name()
{
    return "smargiassi-madden";
}

}
