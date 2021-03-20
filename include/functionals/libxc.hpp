// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "xc.h" // libxc header

#include "functional.hpp"

namespace profess
{

class LibXC : public Functional
{
public:
    LibXC(Box box, std::vector<int> xc_func_ids);
    ~LibXC();
    std::string name() override;
    void set_box(std::array<std::array<double,3>,3> box) override;
    double energy(Double3D den) override;
    std::tuple<double, Double3D> energy_potential(Double3D den) override;
    std::array<std::array<double,3>,3> stress(Double3D den) override;
private:
    std::vector<xc_func_type> _xc_func;
    bool _xc_all_lda;
    void _init(std::vector<int> xc_func_ids);
};

}
