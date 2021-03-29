// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "perdew_burke_ernzerhof.hpp"

namespace profess
{

PerdewBurkeErnzerhof::PerdewBurkeErnzerhof(Box box)
    : LibXC(box, {101, 130})
{
}

std::string PerdewBurkeErnzerhof::name()
{
    return "pbe";
}

}
