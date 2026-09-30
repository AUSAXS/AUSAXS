// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/NeutronFormFactor.h>
#include <utility/Exceptions.h>

#include <array>
#include <numbers>
#include <string>

using namespace ausaxs;
using namespace ausaxs::form_factor;
using neutron::FormFactor;
using constants::atom_t;

namespace {
    // Coherent scattering lengths in fm.
    // Sears, Neutron News 3(3), 26-37 (1992), https://doi.org/10.1080/10448639208218770
    namespace b {
        constexpr double H  = -3.7390;
        constexpr double D  =  6.671;
        constexpr double C  =  6.6460;
        constexpr double N  =  9.36;
        constexpr double O  =  5.803;
        constexpr double S  =  2.847;
        constexpr double Ar =  1.909;
    }

    constexpr double scattering_length(atom_t element) {
        switch (element) {
            case atom_t::C:  return b::C;
            case atom_t::N:  return b::N;
            case atom_t::O:  return b::O;
            case atom_t::S:  return b::S;
            case atom_t::Ar: return b::Ar;
            default: throw except::invalid_argument("form_factor::neutron::scattering_length: No scattering length for element (enum " + std::to_string(static_cast<int>(element)) + ")");
        }
    }

    // X-H distances in Å between the nuclear positions, which are longer than the X-ray distances since the hydrogen electron is displaced into the bond.
    // C-H, N-H, and O-H are the neutron-normalized values of Allen, Acta Cryst. B42, 515-522 (1986).
    // S-H is the neutron value of Allen & Bruno, Acta Cryst. B66, 380-386 (2010).
    constexpr double bond_length(atom_t element) {
        switch (element) {
            case atom_t::C: return 1.083;
            case atom_t::N: return 1.009;
            case atom_t::O: return 0.983;
            case atom_t::S: return 1.338;
            default: throw except::invalid_argument("form_factor::neutron::bond_length: No X-H bond length for element (enum " + std::to_string(static_cast<int>(element)) + ")");
        }
    }

    // H-X-H angles in degrees of the groups with more than one hydrogen.
    // NH2 only occurs in the planar amide and guanidinium groups of proteins.
    constexpr double tetrahedral = 109.4712;
    constexpr double trigonal = 120;
    constexpr double hh_angle(form_factor_t type) {
        switch (type) {
            case form_factor_t::CH2: return tetrahedral;
            case form_factor_t::CH3: return tetrahedral;
            case form_factor_t::NH2: return trigonal;
            case form_factor_t::NH3: return tetrahedral;
            default: throw except::invalid_argument("form_factor::neutron::hh_angle: No H-X-H angle for form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
        }
    }

    constexpr auto make_table(double b_hydrogen) {
        std::array<FormFactor, total_ff_count> table{};
        for (int i = start_index_for_explicit_exv(); i < total_ff_count; ++i) {
            const auto& info = form_factor::detail::ff_info_table[i];
            if (info.element == atom_t::H) {
                table[i] = FormFactor(b_hydrogen, 0, 0, 0, 0);
                continue;
            }
            double d_xh = info.hydrogens == 0 ? 0 : bond_length(info.element);
            double d_hh = info.hydrogens < 2 ? 0 : 2*d_xh*constexpr_math::sin(hh_angle(info.type)*std::numbers::pi/360);
            table[i] = FormFactor(scattering_length(info.element), b_hydrogen, info.hydrogens, d_xh, d_hh);
        }
        return table;
    }

    constexpr auto protonated_table = make_table(b::H);
    constexpr auto deuterated_table = make_table(b::D);

    const FormFactor& get(const std::array<FormFactor, total_ff_count>& table, form_factor_t type) {
        if (type == form_factor_t::EXCLUDED_VOLUME) {
            throw except::runtime_error("form_factor::neutron::get: The excluded volume form factor is not defined for neutrons.");
        }
        if (!form_factor::detail::is_tabulated(type)) {
            throw except::runtime_error("form_factor::neutron::get: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
        }
        return table[static_cast<int>(type)];
    }
}

const FormFactor& neutron::protonated::get(form_factor_t type) {
    return ::get(protonated_table, type);
}

const FormFactor& neutron::deuterated::get(form_factor_t type) {
    return ::get(deuterated_table, type);
}
