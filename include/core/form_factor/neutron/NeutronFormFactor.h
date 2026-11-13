// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsProperties.h>
#include <form_factor/FormFactorInfo.h>
#include <form_factor/GroupFormFactor.h>
#include <form_factor/neutron/NeutronFormFactorTable.h>
#include <utility/Exceptions.h>

#include <array>
#include <string>
#include <utility>

namespace ausaxs::form_factor::neutron {
    /**
     * @brief The neutron form factor of a single nucleus: its coherent scattering length, which does not depend on q.
     */
    class AtomicFormFactor {
        public:
            constexpr AtomicFormFactor(double scattering_length) : b(scattering_length) {}

            /**
             * @brief Evaluate the form factor at a given q value.
             */
            constexpr double evaluate(double /*q*/) const {return b;}

            /**
             * @brief Evaluate the form factor at q = 0.
             */
            constexpr double I0() const {return b;}

        private:
            double b;
    };

    /**
     * @brief The neutron form factor of a form factor type: a heavy atom with its bound hydrogens.
     */
    using FormFactor = GroupFormFactor<AtomicFormFactor>;

    namespace detail {
        // the sample is assumed to be protonated and in 100% H2O, so every hydrogen is protium
        constexpr AtomicFormFactor hydrogen(scattering_length::H);

        constexpr auto table = [] <std::size_t... I> (std::index_sequence<I...>) {
            return std::array<FormFactor, total_ff_count>{FormFactor(
                AtomicFormFactor(form_factor::detail::ff_info_table[I].neutron_scattering_length),
                hydrogen,
                form_factor::detail::ff_info_table[I].geometry
            )...};
        }(std::make_index_sequence<total_ff_count>{});
    }

    /**
     * @brief Get the neutron form factor of a form factor type, as described by form_factor::detail::ff_info_table.
     *        All hydrogens, implicit and explicit, are protium.
     *        The excluded volume is not defined for neutrons, and requesting it throws.
     */
    constexpr const FormFactor& get(form_factor_t type) {
        if (type == form_factor_t::EXCLUDED_VOLUME) {
            throw ausaxs::except::runtime_error("form_factor::neutron::get: The excluded volume form factor is not defined for neutrons.");
        }
        if (!form_factor::detail::is_tabulated(type)) {
            throw ausaxs::except::runtime_error("form_factor::neutron::get: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
        }
        return detail::table[static_cast<int>(type)];
    }

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
     */
    constexpr double solvent_density(double d2o_fraction) {
        if (d2o_fraction < 0 || 1 < d2o_fraction) {
            throw ausaxs::except::invalid_argument("form_factor::neutron::solvent_density: The D2O fraction must be between 0 and 1, but was " + std::to_string(d2o_fraction) + ".");
        }
        return (1-d2o_fraction)*detail::density_H2O + d2o_fraction*detail::density_D2O;
    }
}
