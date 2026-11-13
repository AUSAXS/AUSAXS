// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

// Neutron scattering lengths. The geometry of the united atomic groups is shared by all probes, see form_factor/GroupGeometry.h.
/**
 * @brief Coherent scattering lengths in fm, as distributed by periodictable 2.0.2 (also used by SasView).
 *        Natural elements are from Rauch & Waschkowski, ILL Neutron Data Booklet, 2nd ed. (2003),
 *        with newer measurements for 1H, 2H, C, O, Sn, Pb, Bi from Snow et al., Phys. Rev. D 101, 062004 (2020), https://doi.org/10.1103/PhysRevD.101.062004
 *        and for 4He from Haun et al., Phys. Rev. Lett. 124, 012501 (2020), https://doi.org/10.1103/PhysRevLett.124.012501 (natural He is its isotope-weighted average with 3He).
 *        H is protium, not natural hydrogen, since every hydrogen is treated as protium.
 */
namespace ausaxs::form_factor::neutron::scattering_length {
    constexpr double H  = -3.7395;
    constexpr double D  =  6.6681;
    constexpr double He =  3.0985;
    constexpr double Li =    -1.9;
    constexpr double Be =    7.79;
    constexpr double B  =     5.3;
    constexpr double C  =  6.6472;
    constexpr double N  =    9.36;
    constexpr double O  =  5.8037;
    constexpr double F  =   5.654;
    constexpr double Ne =   4.566;
    constexpr double Na =    3.63;
    constexpr double Mg =   5.375;
    constexpr double Al =   3.449;
    constexpr double Si = 4.15071;
    constexpr double P  =    5.13;
    constexpr double S  =   2.847;
    constexpr double Cl =  9.5792;
    constexpr double Ar =   1.909;
    constexpr double K  =    3.67;
    constexpr double Ca =     4.7;
    constexpr double Sc =    12.1;
    constexpr double Ti =   -3.37;
    constexpr double V  =  -0.443;
    constexpr double Cr =   3.635;
    constexpr double Mn =   -3.75;
    constexpr double Fe =    9.45;
    constexpr double Co =    2.49;
    constexpr double Ni =    10.3;
    constexpr double Cu =   7.718;
    constexpr double Zn =    5.68;
    constexpr double Ga =   7.288;
    constexpr double Ge =   8.185;
    constexpr double As =    6.58;
    constexpr double Se =    7.97;
    constexpr double Br =    6.79;
    constexpr double Kr =    7.81;
    constexpr double Rb =    7.08;
    constexpr double Sr =    7.02;
    constexpr double Y  =    7.75;
    constexpr double Zr =    7.16;
    constexpr double Nb =   7.054;
    constexpr double Mo =   6.715;
    constexpr double Tc =     6.8;
    constexpr double Ru =    7.02;
    constexpr double Rh =     5.9;
    constexpr double Pd =    5.91;
    constexpr double Ag =   5.922;
    constexpr double Cd =    4.83;
    constexpr double In =   4.065;
    constexpr double Sn =  6.2239;
    constexpr double Sb =    5.57;
    constexpr double Te =    5.68;
    constexpr double I  =    5.28;
    constexpr double Xe =    4.69;
    constexpr double Cs =    5.42;
    constexpr double Ba =    5.07;
    constexpr double Hf =    7.77;
    constexpr double Ta =    6.91;
    constexpr double W  =   4.755;
    constexpr double Re =     9.2;
    constexpr double Os =    10.7;
    constexpr double Ir =    10.6;
    constexpr double Pt =     9.6;
    constexpr double Au =     7.9;
    constexpr double Tl =   8.776;
    constexpr double Pb =  9.4024;
    constexpr double Bi =  8.5242;
    constexpr double other = Ar; // default scattering length for unknown atoms

    // the excluded volume has no scattering length
    constexpr double excluded_volume = 0;
}
