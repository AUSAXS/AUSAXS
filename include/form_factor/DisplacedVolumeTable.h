#pragma once

#include <constants/SI.h>
#include <constants/vdwTable.h>
#include <math/ConstMath.h>

#include <cmath>

#define TRAUBE_FF true
#define PONTIUS_FF true
namespace constants::displaced_volume {
    namespace {
        constexpr double volume(double radius) {
            return 4*M_PI/3*math::pow(radius, 3);
        }
    }

    #if TRAUBE_FF
        // from original CRYSOL paper: https://doi.org/10.1107/S0021889895007047
        constexpr double CH  = 0.02159*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double CH2 = 0.02674*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double CH3 = 0.03189*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double NH  = 0.00764*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double NH2 = 0.01279*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double NH3 = 0.01794*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double OH  = 0.01428*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double SH  = 0.02510*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);

        constexpr double H = 0.00515*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double C = 0.01644*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double N = 0.00249*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double O = 0.00913*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double S = 0.01986*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);

    #elif PONTIUS_FF
        // from Pontius et al, 1996: https://doi.org/10.1006/jmbi.1996.0628
        constexpr double CH = 11.8;
        constexpr double CH2 = 20.9;
        constexpr double CH3 = 33.9;
        constexpr double NH = 14.1;
        constexpr double NH2 = 24.6;
        constexpr double NH3 = 23;
        constexpr double OH = 23;
        constexpr double SH = 34.2;

        constexpr double H = 0.00515*math::pow(constants::SI::length::nm/constants::SI::length::A, 3);
        constexpr double C = 8.4;
        constexpr double N = 8.8;
        constexpr double O = 22.3;
        constexpr double S = 25;

    #else 
        constexpr double CH = volume(constants::radius::vdw::C) + volume(constants::radius::vdw::H);
        constexpr double CH2 = volume(constants::radius::vdw::C) + 2*volume(constants::radius::vdw::H);
        constexpr double CH3 = volume(constants::radius::vdw::C) + 3*volume(constants::radius::vdw::H);
        constexpr double NH = volume(constants::radius::vdw::N) + volume(constants::radius::vdw::H);
        constexpr double NH2 = volume(constants::radius::vdw::N) + 2*volume(constants::radius::vdw::H);
        constexpr double NH3 = volume(constants::radius::vdw::N) + 3*volume(constants::radius::vdw::H);
        constexpr double OH = volume(constants::radius::vdw::O) + volume(constants::radius::vdw::H);
        constexpr double SH = volume(constants::radius::vdw::S) + volume(constants::radius::vdw::H);

        constexpr double H = volume(constants::radius::vdw::H);
        constexpr double C = volume(constants::radius::vdw::C);
        constexpr double N = volume(constants::radius::vdw::N);
        constexpr double O = volume(constants::radius::vdw::O);
        constexpr double S = volume(constants::radius::vdw::S);
    #endif
    constexpr double avg_vol = 50./3;
    constexpr double OH2 = 2.98*math::pow(10, -23)*math::pow(constants::SI::length::cm/constants::SI::length::A, 3);
    constexpr double Ar = volume(constants::radius::vdw::Ar);
}