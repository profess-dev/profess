// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "libxc.hpp"

namespace profess
{

class PerdewBurkeErnzerhof : public LibXC
{
public:
    PerdewBurkeErnzerhof(Box box);
    std::string name();
};

}
