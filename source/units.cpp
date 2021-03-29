// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "units.hpp"

#include <stdexcept>

namespace profess
{

namespace units
{

double convert_length(double length, std::string in, std::string out)
{
    if (in == "b") {
        if (out == "b") return length;
        if (out == "m") return length * _to_meters;
        if (out == "a") return length * _to_angstroms;
    } else if (in == "m") {
        if (out == "b") return length / _to_meters;
        if (out == "m") return length;
        if (out == "a") return length / _to_meters * _to_angstroms;
    } else if (in == "a") {
        if (out == "b") return length / _to_angstroms;
        if (out == "m") return length / _to_angstroms * _to_meters;
        if (out == "a") return length;
    }

    // should have returned by now
    throw std::runtime_error("units error");
}

std::array<std::array<double,3>,3> convert_length(
        std::array<std::array<double,3>,3> lengths,
        std::string in,
        std::string out)
{
    auto new_lengths = lengths;
    for (auto& array : new_lengths) {
        for (auto& element : array) {
            element = convert_length(element, in, out);
        }
    }
    return new_lengths;
}

std::vector<std::array<double,3>> convert_length(
        std::vector<std::array<double,3>> lengths,
        std::string in,
        std::string out)
{
    auto new_lengths = lengths;
    for (auto& array : new_lengths) {
        for (auto& element : array) {
            element = convert_length(element, in, out);
        }
    }
    return new_lengths;
}

double convert_volume(double volume, std::string in, std::string out)
{
    if (in == "b3") {
        if (out == "b3") return volume;
        if (out == "m3") return volume * _to_meters3;
        if (out == "a3") return volume * _to_angstroms3;
    } else if (in == "m3") {
        if (out == "b3") return volume / _to_meters3;
        if (out == "m3") return volume;
        if (out == "a3") return volume / _to_meters3 * _to_angstroms3;
    } else if (in == "a3") {
        if (out == "b3") return volume / _to_angstroms3;
        if (out == "m3") return volume / _to_angstroms3 * _to_meters3;
        if (out == "a3") return volume;
    }

    // should have returned by now
    throw std::runtime_error(std::string("in convert_volume, can't convert ")
                                + in + std::string(" to ") + out);
}

double convert_energy(double energy, std::string in, std::string out)
{
    if (in == "h") {
        if (out == "h") return energy;
        if (out == "j") return energy * _to_joules;
        if (out == "ev") return energy * _to_electron_volts;
    } else if (in == "j") {
        if (out == "h") return energy / _to_joules;
        if (out == "j") return energy;
        if (out == "ev") return energy / _to_joules * _to_electron_volts;
    } else if (in == "ev") {
        if (out == "h") return energy / _to_electron_volts;
        if (out == "j") return energy / _to_electron_volts * _to_joules;
        if (out == "ev") return energy;
    }

    // should have returned by now
    throw std::runtime_error("units error");
}

double convert_pressure(double pressure, std::string in, std::string out)
{
    if (in == "h/b3") {
        if (out == "h/b3") return pressure;
        if (out == "gpa") return pressure * _to_gigapascals;
        if (out == "ev/a3") return pressure * _to_electron_volts_per_angstrom3;
    } else if (in == "gpa") {
        if (out == "h/b3") return pressure / _to_gigapascals;
        if (out == "gpa") return pressure;
        if (out == "ev/a3") return pressure / _to_gigapascals
                                        * _to_electron_volts_per_angstrom3;
    } else if (in == "ev/a3") {
        if (out == "h/b3") return pressure / _to_electron_volts_per_angstrom3;
        if (out == "gpa") return pressure / _to_electron_volts_per_angstrom3
                                   * _to_gigapascals;
        if (out == "ev/a3") return pressure;
    }

    // should have returned by now
    throw std::runtime_error("units error");
}

std::array<std::array<double,3>,3> convert_pressure(
        std::array<std::array<double,3>,3> pressures,
        std::string in,
        std::string out)
{
    auto new_pressures = pressures;
    for (auto& array : new_pressures) {
        for (auto& element : array) {
            element = convert_pressure(element, in, out);
        }
    }
    return new_pressures;
}

}

}
