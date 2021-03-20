// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include "functional.hpp"

namespace profess
{

class KineticGGA : public Functional
{
public:
    KineticGGA(
        Box box,
        std::function<std::tuple<Double3D,Double3D,Double3D>
                      (Double3D,Double3D)> integrand_and_derivatives,
        std::function<std::tuple<Double3D,Double3D>
                      (Double3D,Double3D)> second_derivatives);
    std::string name() override;
    void set_box(std::array<std::array<double,3>,3> box) override;
    double energy(Double3D den) override;
    std::tuple<double,Double3D> energy_potential(Double3D den) override;
    std::array<std::array<double,3>,3> stress(Double3D den) override;
private:
    std::function<std::tuple<Double3D,Double3D,Double3D>
                  (Double3D,Double3D)> _integrand_and_derivatives;
    std::function<std::tuple<Double3D,Double3D>
                  (Double3D,Double3D)> _second_derivatives;
};

}
