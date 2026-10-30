// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

// Neutron scattering lengths. The geometry of the united atomic groups is shared by all probes, see form_factor/GroupGeometry.h.
/**
 * @brief Coherent scattering lengths in fm.
 *        Sears, Neutron News 3(3), 26-37 (1992), https://doi.org/10.1080/10448639208218770
 */
namespace ausaxs::form_factor::neutron::scattering_length {
    constexpr double H  = -3.7390;
    constexpr double D  =  6.671;
    constexpr double C  =  6.6460;
    constexpr double N  =  9.36;
    constexpr double O  =  5.803;
    constexpr double F  =  5.654;
    constexpr double Na =  3.63;
    constexpr double Mg =  5.375;
    constexpr double P  =  5.13;
    constexpr double S  =  2.847;
    constexpr double Cl =  9.5770;
    constexpr double Ar =  1.909;
    constexpr double K  =  3.67;
    constexpr double Ca =  4.70;
    constexpr double Mn = -3.73;
    constexpr double Fe =  9.45;
    constexpr double Co =  2.49;
    constexpr double Ni = 10.3;
    constexpr double Cu =  7.718;
    constexpr double Zn =  5.680;
    constexpr double Se =  7.970;
    constexpr double Br =  6.795;
    constexpr double I  =  5.28;
    constexpr double other = Ar; // default scattering length for unknown atoms

    // the excluded volume has no scattering length
    constexpr double excluded_volume = 0;
}
