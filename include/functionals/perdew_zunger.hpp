// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "libxc.hpp"

namespace profess
{

class PerdewZunger : public LibXC
{
public:
    PerdewZunger(Box box);
    std::string name();
};

}
