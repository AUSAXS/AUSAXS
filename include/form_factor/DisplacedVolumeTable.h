#pragma once

#include <form_factor/FormFactorType.h>
#include <constants/Constants.h>
#include <constants/vdwTable.h>
#include <math/ConstMath.h>

#include <cmath>

namespace constants::DisplacedVolume {
    namespace {
        constexpr double volume(double radius) {
            return 4*M_PI/3*math::pow(radius, 3);
        }
    }

    // from original CRYSOL paper: https://doi.org/10.1107/S0021889895007047
    constexpr double CH =  0.02159*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double CH2 = 0.02674*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double CH3 = 0.03189*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double NH =  0.00764*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double NH2 = 0.01279*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double NH3 = 0.01794*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double OH =  0.01428*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    // constexpr double SH =  0.02510*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double SH = 8.35334/0.334;

    constexpr double H = 0.00515*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double C = 0.01644*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double N = 0.00249*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double O = 0.00913*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
    constexpr double S = 0.01986*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);

    // constexpr double H =  volume(constants::radius::vdw::H);
    // constexpr double C =  volume(constants::radius::vdw::C);
    // constexpr double N =  volume(constants::radius::vdw::N);
    // constexpr double O =  volume(constants::radius::vdw::O);
    // constexpr double S =  volume(constants::radius::vdw::S);
    constexpr double Ar = volume(constants::radius::vdw::Ar);
}