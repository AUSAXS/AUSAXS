// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/tables/FormFactorTableHelper.h>

// Neutron scattering lengths in units of fm. 
namespace ausaxs::constants::form_factor::neutron {

    // https://www.nist.gov/ncnr/neutron-scattering-lengths-list
    constexpr double H = -3.7390;
    constexpr double C = 6.6460;
    constexpr double N = 9.36;
    constexpr double O = 5.803;
    constexpr double S = 2.847;
    constexpr double Ar = 1.909;
    constexpr double other = Ar;

    // unified groups, which are trivially calculated as the sum of their parts due to the q independence
    constexpr double CH = C + H;
    constexpr double CH2 = C + 2*H;
    constexpr double CH3 = C + 3*H;
    constexpr double NH = N + H;
    constexpr double NH2 = N + 2*H;
    constexpr double NH3 = N + 3*H;
    constexpr double OH = O + H;
    constexpr double SH = S + H;
}