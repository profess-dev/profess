// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "wang_govind_carter_1999_i.hpp"

namespace profess
{

WangGovindCarter1999I::WangGovindCarter1999I(
        Box box,
        std::array<size_t,3> grid_shape,
        double den0,
        bool den0_immutable)
    : KineticClassA(
        box,
        grid_shape,
        (5.0+std::sqrt(5.0))/6.0,
        (5.0-std::sqrt(5.0))/6.0,
        [](double x){return 1.0+x;},
        [](double x){return 1.0;},
        den0,
        den0_immutable)
{
}

std::string WangGovindCarter1999I::name()
{
    return "wang-govind-carter-1999-i";
}

}
