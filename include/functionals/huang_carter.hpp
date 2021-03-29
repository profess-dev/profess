// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once
#include "libKEDF.h"
#include "functional.hpp"

namespace profess
{

class HuangCarter : public Functional
{
public:
    HuangCarter(
        Box box,
        std::array<size_t,3> grid_shape,
        double den0,
        double lambda=0.01,
        double alpha=2.0166666666666666,
        double beta=0.65,
        double refRatio=1.15);
    std::string name() override;
    void set_box(std::array<std::array<double,3>,3> box) override;
    double energy(Double3D den) override;
    std::tuple<double,Double3D> energy_potential(Double3D den) override;
    std::array<std::array<double,3>,3> stress(Double3D den) override;
private:
    libkedf_data* _kedf_data;
    double _den0;
    const double _lambda;
    const double _alpha;
    const double _beta;
    const double _refRatio;
};

}
