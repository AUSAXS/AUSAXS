// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsProperties.h>
#include <form_factor/neutron/NeutronFormFactorTable.h>
#include <utility/Exceptions.h>

#include <string>

/**
 * @brief Conversions between the neutron scattering length density of an H2O/D2O solvent and its D2O fraction.
 *        settings::scattering::neutron::solvent_density is the single source of truth for the solvent.
 *        These helpers only translate between it and the H2O/D2O mixture that realises it.
 */
namespace ausaxs::form_factor::neutron {
    namespace detail {
        // Number of molecules per Å^3 of a liquid with the given mass density (Da/Å^3) and molecular mass (Da).
        constexpr double number_density(double mass_density, double molecular_mass) {
            return mass_density/molecular_mass;
        }

        constexpr double density_H2O = (2*scattering_length::H + scattering_length::O)*number_density(
            constants::mass::density::water, 2*constants::mass::get_mass(constants::atom_t::H) + constants::mass::get_mass(constants::atom_t::O)
        );
        constexpr double density_D2O = (2*scattering_length::D + scattering_length::O)*number_density(
            constants::mass::density::D2O, 2*constants::mass::get_mass(constants::atom_t::D) + constants::mass::get_mass(constants::atom_t::O)
        );
    }

    /**
     * @brief The coherent scattering length density of an H2O/D2O mixture in fm/Å^3.
     *
     * @param d2o_fraction The D2O fraction of the mixture, between 0 and 1.
     */
    constexpr double solvent_density(double d2o_fraction) {
        if (d2o_fraction < 0 || 1 < d2o_fraction) {
            throw ausaxs::except::invalid_argument("form_factor::neutron::solvent_density: The D2O fraction must be between 0 and 1, but was " + std::to_string(d2o_fraction) + ".");
        }
        return (1-d2o_fraction)*detail::density_H2O + d2o_fraction*detail::density_D2O;
    }

    /**
     * @brief The D2O fraction of the H2O/D2O mixture with the given coherent scattering length density in fm/Å^3. This is the inverse of solvent_density.
     *        A density outside the range spanned by pure H2O and pure D2O cannot be realised by water alone.
     *        It is extrapolated linearly, so the result then lies outside [0, 1].
     */
    constexpr double d2o_fraction(double solvent_density) {
        return (solvent_density - detail::density_H2O)/(detail::density_D2O - detail::density_H2O);
    }
}
