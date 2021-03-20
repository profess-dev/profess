// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include <functional>
#include <string>
#include <vector>

namespace profess
{

class Ions
{
public:

    size_t count() const;
    size_t count_types() const;
    size_t count_of_type(size_t i) const;
    double count_charge() const;

    std::vector<std::array<double,3>>& xyz_coords();
    const std::vector<std::array<double,3>>& xyz_coords() const;

    const std::vector<size_t>& type_ids() const;

    std::vector<std::array<double,3>> xyz_coords_of_type(size_t i) const;

    // type info
    const std::vector<double>& charges() const;
    const std::vector<std::function<double(double)>>& ft_potentials() const;
    std::vector<std::function<double(double)>>
            ft_potential_derivatives();

    Ions& add_ion_type_coulomb(double z);
    Ions& add_ion_type_recpot(std::string filename);
    Ions& add_ion_type_harmonic_compactified(double w, double r1, double r2);
    Ions& add_ion_type_generic(
            std::function<double(double)> ft_potential,
            std::function<double(double)> ft_potential_derivative,
            double charge=0.0);

    Ions& set_ions(
            std::vector<std::array<double,3>> xyz_coords,
            std::vector<size_t> type_ids);

private:

    std::vector<std::array<double,3>> _xyz_coords;
    std::vector<size_t> _type_ids;

    // arrays for ion types
    std::vector<double> _charges;
    std::vector<std::function<double(double)>> _ft_potentials;
    std::vector<std::function<double(double)>> _ft_potential_derivatives;

};

}
