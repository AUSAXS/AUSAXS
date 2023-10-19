#pragma once

#include <form_factor/FormFactorType.h>
#include <constants/Constants.h>
#include <constants/vdwTable.h>

#include <cmath>

namespace constants::DisplacedVolume {
    namespace {
        constexpr double volume(double radius) {
            return 4*M_PI/3*std::pow(radius, 3);
        }
    }

    // from original CRYSOL paper: https://doi.org/10.1107/S0021889895007047
    constexpr double CH =  0.02159*std::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double CH2 = 0.02674*std::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double CH3 = 0.03189*std::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double NH =  0.00764*std::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double NH2 = 0.01279*std::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double NH3 = 0.01794*std::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double OH =  0.01428*std::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double SH =  0.02510*std::pow(constants::SI::length::nm/constants::SI::length::A, 3);

    constexpr double H =  volume(constants::radius::vdw::H);
    constexpr double C =  volume(constants::radius::vdw::C);
    constexpr double N =  volume(constants::radius::vdw::N);
    constexpr double O =  volume(constants::radius::vdw::O);
    constexpr double S =  volume(constants::radius::vdw::S);
    constexpr double Ar = volume(constants::radius::vdw::Ar);
}