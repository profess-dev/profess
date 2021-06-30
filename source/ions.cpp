// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "ions.hpp"
#include "tools.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

namespace profess
{

size_t Ions::count() const
{
    return _xyz_coords.size();
}

size_t Ions::count_types() const
{
    return _ft_potentials.size();
}

size_t Ions::count_of_type(size_t i) const
{
    size_t count = 0;
    for (size_t ii=0; ii<_type_ids.size(); ++ii)
        if (_type_ids[ii] == i)
            count++;
    return count;
}

double Ions::count_charge() const
{
    double charge = 0.0;
    for (size_t i=0; i<count_types(); ++i) {
        charge += count_of_type(i) * _charges[i];
    }
    return charge;
}

std::vector<std::array<double,3>>& Ions::xyz_coords()
{
    return _xyz_coords;
}

const std::vector<std::array<double,3>>& Ions::xyz_coords() const
{
    return _xyz_coords;
}

const std::vector<size_t>& Ions::type_ids() const
{
    return _type_ids;
}


std::vector<std::array<double,3>> Ions::xyz_coords_of_type(size_t i) const
{
    std::vector<std::array<double,3>> xyz_coords(count_of_type(i));
    size_t j = 0;
    for (size_t ii=0; ii<count(); ++ii) {
        if (_type_ids[ii] == i) {
            xyz_coords[j] = _xyz_coords[ii];
            j += 1;
        }
    }
    return xyz_coords;
}


const std::vector<double>& Ions::charges() const
{
    return _charges;
}

const std::vector<std::function<double(double)>>& Ions::ft_potentials() const
{
    return _ft_potentials;
}

std::vector<std::function<double(double)>>
        Ions::ft_potential_derivatives()
{
    return _ft_potential_derivatives;
}

Ions& Ions::add_ion_type_coulomb(double z, double cutoff)
{
    _charges.emplace_back(z);
    // no cutoff
    if (cutoff < 0.0) {
        _ft_potentials.emplace_back(
            [z](double k) {
                if (k<1e-8) {
                    return 0.0;
                } else {
                    return -4.0*M_PI*z/(k*k);
                }
            }
        );
        _ft_potential_derivatives.emplace_back(
            [z](double k) {
                if (k<1e-8) {
                    return 0.0;
                } else {
                    return 8.0*M_PI*z/(k*k*k);
                }
            }
        );
    // with cutoff
    } else {
        _ft_potentials.emplace_back(
            [z,cutoff](double k) {
                if (k<1e-8) {
                    return -2.0*M_PI*z*cutoff*cutoff;
                } else {
                    return -4.0*M_PI*z/(k*k)*(1.0-std::cos(k*cutoff));
                }
            }
        );
        _ft_potential_derivatives.emplace_back(
            [z,cutoff](double k) {
                if (k<1e-8) {
                    return 0.0;
                } else {
                    return 8.0*M_PI*z/(k*k*k)*(
                        1.0-std::cos(k*cutoff)-0.5*k*cutoff*std::sin(k*cutoff));
                }
            }
        );
    }
    return *this;
}

Ions& Ions::add_ion_type_recpot(std::string filename)
{
    const double _bohr = 0.529177208607388;
    const double _hartree_to_ev = 27.2113834279111; 
    // open file
    std::ifstream f(filename);
    // count lines
    std::string line;
    size_t num_lines = 0;
    while (std::getline(f, line)) num_lines++;
    f.clear();
    f.seekg(0);
    // extract comment block (ignored)
    line = "";
    size_t num_header_lines = 0;
    while (line.compare(0, 11, "END COMMENT", 0, 11)) {
        std::getline(f, line);
        ++num_header_lines;
    }
    // extract version (ignored)
    std::getline(f, line);
    // extract k_max in [1/angstrom] and convert to [1/bohr]
    double k_max;
    f >> k_max;
    k_max *= _bohr;
    // extract potential and convert to atomic units
    size_t size = 3*(num_lines - num_header_lines - 2 - 1);
    std::vector<double> ft_pot_data(size);
    double to_au = 1.0 / (_bohr * _bohr * _bohr * _hartree_to_ev);
    for (size_t i=0; i<size; ++i) {
        f >> ft_pot_data[i];
        ft_pot_data[i] *= to_au;
    }
    // check file terminator
    double value;
    f >> value;
    if (value != 1000) {
        std::runtime_error("recpot file terminated improperly.\n");
    }
    // compute smallest nonzero wavevector
    double k_inc = k_max / (ft_pot_data.size() - 1);
    // extract charge from coulombic part of potential
    double z = ft_pot_data[1] - ft_pot_data[0]; // get coulombic part at [1]
    z *= k_inc * k_inc / (-4.0 * M_PI); // extract charge
    z = std::round(z); // round to nearest integer
    // remove coulombic part from ft_pot_data
    for (size_t a=1; a<size; ++a) {
        ft_pot_data[a] += 4.0 * M_PI * z / ( a*k_inc * a*k_inc);
    }
    // set ion charge
    _charges.emplace_back(z);
    // set the ion potential function
    _ft_potentials.emplace_back(
        [ft_pot_data, k_inc, k_max, z](double k) {
            if (k > k_max) {
                return 0.0;
            } else {
                double pot = interpolate_cubic_lagrange(ft_pot_data, k, 0.0, k_inc);
                if (k < 1e-12) {
                    return pot; // add coulombic part only for k>0
                } else {
                    return pot - 4.0*M_PI*z/(k*k);
                }
            }
        }
    );
    // set the ion potential derivative
    _ft_potential_derivatives.emplace_back(
        [ft_pot_data, k_inc, k_max, z](double k) {
            if (k > k_max) {
                return 0.0;
            } else {
                double pot = interpolate_derivative_cubic_lagrange(
                        ft_pot_data, k, 0.0, k_inc);
                if (k < 1e-12) {
                    return pot; // add coulombic part only for k>0
                } else {
                    return pot + 8.0*M_PI*z/(k*k*k);
                }
            }
        }
    );
    return *this;
}

Ions& Ions::add_ion_type_harmonic_compactified(double w, double r1, double r2)
{
    // originates from 1/2*r^2 attenuated by smoothstep function
    _charges.push_back(0.0);
    _ft_potentials.push_back(
        [w,r1,r2](double k) {
            if (k<1e-8) {
                const double r17 = r1*r1*r1*r1*r1*r1*r1;
                const double r27 = r2*r2*r2*r2*r2*r2*r2;
                return M_PI*w*w/(70*(r1-r2)*(r1-r2)*(r1-r2))
                        *(3*r1*r17-4*r17*r2+4*r1*r27-3*r2*r27);
            } else {
                const double k2 = k*k;
                const double k4 = k2*k2;
                return 12*M_PI*w*w/(k4*k4*(r1-r2)*(r1-r2)*(r1-r2))
                    *(  ( 240+k4*r1*r1*r1*(r1-r2)
                                 +12*k2*r1*(-5*r1+3*r2))*std::cos(k*r1)
                       +(-240+k4*(r1-r2)*r2*r2*r2
                                 +12*k2*r2*(-3*r1+5*r2))*std::cos(k*r2)
                       +( 180*r1-11*k2*r1*r1*r1
                                 -60*r2+9*k2*r1*r1*r2)*k*std::sin(k*r1) 
                       +(-180*r2+11*k2*r2*r2*r2
                                 +60*r1-9*k2*r1*r2*r2)*k*std::sin(k*r2) );
            }
        }
    );
    _ft_potential_derivatives.push_back(
        [w,r1,r2] (double k) {
            if (k<1e-8) {
                return 0.0;
            } else {
                const double k2 = k*k;
                const double k4 = k2*k2;
                return -12*M_PI*w*w/(k4*k4*k*(r1-r2)*(r1-r2)*(r1-r2))
                    *(  ( 1920+k4*r1*r1*r1*(15*r1-13*r2)
                                 +12*k2*r1*(-45*r1+23*r2))*std::cos(k*r1)
                       +(-1920+k4*(13*r1-15*r2)*r2*r2*r2
                                 +12*k2*r2*(-23*r1+45*r2))*std::cos(k*r2)
                       +( 1500*r1-115*k2*r1*r1*r1+k4*r1*r1*r1*r1*(r1-r2)
                                 -420*r2+81*k2*r1*r1*r2)*k*std::sin(k*r1) 
                       +(-1500*r2+115*k2*r2*r2*r2+k4*r2*r2*r2*r2*(r1-r2)
                                 +420*r1-81*k2*r1*r2*r2)*k*std::sin(k*r2) );
            }
        }
    );
    return *this;
}

Ions& Ions::add_ion_type_generic(
        std::function<double(double)> ft_potential,
        std::function<double(double)> ft_potential_derivative,
        double charge)
{
    _charges.push_back(charge);
    _ft_potentials.push_back(ft_potential);
    _ft_potential_derivatives.push_back(ft_potential_derivative);
    return *this;
}

Ions& Ions::set_ions(
        std::vector<std::array<double,3>> xyz_coords,
        std::vector<size_t> type_ids)
{
    // check that vectors holding xyz_coords and type_ids have equal lengths
    if (xyz_coords.size() != type_ids.size()) {
        throw std::runtime_error("xyz_coords and type_ids "
                                 "should be the same size.");
    }
    // check that type ids begin with 0
    size_t min_type_id = *std::min_element(type_ids.begin(), type_ids.end());
    if (min_type_id != 0) {
        throw std::runtime_error("type_ids should begin with zero.");
    }
    // set xyz_coords and type_ids
    _xyz_coords = xyz_coords;
    _type_ids = type_ids;
    return *this;
}

}
