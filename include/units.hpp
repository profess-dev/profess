// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <string>
#include <vector>

namespace profess
{

namespace units
{

// ----- based on 2018 CODATA recommended values -----

// length unit: bohr radius
const double _to_meters = 5.29177210903e-11;
const double _to_angstroms = _to_meters * 1e10;
double convert_length(double length, std::string in, std::string out);
std::array<std::array<double,3>,3> convert_length(
        std::array<std::array<double,3>,3> lengths,
        std::string in,
        std::string out);
std::vector<std::array<double,3>> convert_length(
        std::vector<std::array<double,3>> lengths,
        std::string in,
        std::string out);

// volume unit: bohr^3
const double _to_meters3 = _to_meters * _to_meters * _to_meters;
const double _to_angstroms3 = _to_angstroms * _to_angstroms * _to_angstroms;
double convert_volume(double volume, std::string in, std::string out);

// energy unit: hartree
const double _to_joules = 4.3597447222071e-18;
const double _to_electron_volts = _to_joules / 1.602176634e-19;
double convert_energy(double energy, std::string in, std::string out);

// pressure unit: hartree/bohr^3
const double _to_gigapascals =
        _to_joules / (_to_meters * _to_meters * _to_meters) * 1e-9;
const double _to_electron_volts_per_angstrom3 =
        _to_electron_volts / (_to_angstroms * _to_angstroms * _to_angstroms);
double convert_pressure(double pressure, std::string in, std::string out);
std::array<std::array<double,3>,3> convert_pressure(
        std::array<std::array<double,3>,3> pressures,
        std::string in,
        std::string out);

// force unit: hartree/bohr
const double _to_electron_volts_per_angstrom = _to_electron_volts / _to_angstroms;
double convert_force(double force, std::string in, std::string out);
std::vector<std::array<double,3>> convert_force(
        std::vector<std::array<double,3>> forces,
        std::string in,
        std::string out);

}

}
