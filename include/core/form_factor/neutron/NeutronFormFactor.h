// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsProperties.h>
#include <form_factor/FormFactorInfo.h>
#include <form_factor/neutron/NeutronFormFactorTable.h>
#include <math/ConstexprMath.h>
#include <utility/Exceptions.h>

#include <array>
#include <string>
#include <utility>

namespace ausaxs::form_factor::neutron {
    /**
     * @brief The neutron form factor of a heavy atom X with n equivalent bound hydrogens.
     */
    class FormFactor {
        public:
            /**
             * @brief Initialize a form factor from a tabulated group.
             */
            constexpr FormFactor(const coefficients::Group& group) : group(group) {}

            /**
             * @brief Evaluate the amplitude f(q), which enters every pair of distinct groups as f_i(q)*f_j(q)*sinc(q*r_ij).
             */
            constexpr double evaluate(double q) const {
                return group.b_heavy + group.hydrogens*b_hydrogen*math::sinc(q*group.d_xh);
            }

            /**
             * @brief Evaluate the self-term s(q), which replaces f(q)^2 for a group paired with itself.
             */
            constexpr double evaluate_self(double q) const {
                const auto& [b_heavy, n, d_xh, d_hh] = group;
                return b_heavy*b_heavy + n*b_hydrogen*b_hydrogen
                    + 2*n*b_heavy*b_hydrogen*math::sinc(q*d_xh)
                    + n*(n-1)*b_hydrogen*b_hydrogen*math::sinc(q*d_hh);
            }

            /**
             * @brief Evaluate the amplitude at q = 0.
             */
            constexpr double I0() const {
                return group.b_heavy + group.hydrogens*b_hydrogen;
            }

        private:
            coefficients::Group group;

            // the sample is assumed to be protonated and in 100% H2O, so every hydrogen is protium
            static constexpr double b_hydrogen = scattering_length::H;
    };

    namespace detail {
        constexpr auto table = [] <std::size_t... I> (std::index_sequence<I...>) {
            return std::array<FormFactor, total_ff_count>{FormFactor(form_factor::detail::ff_info_table[I].neutron_coefficients)...};
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
        // Number of molecules per Å^3 of a liquid with the given mass density (SI) and molecular mass (u).
        constexpr double number_density(double mass_density, double molecular_mass) {
            return mass_density*constants::SI::volume::A3/(molecular_mass*constants::SI::mass::u);
        }

        // Both mass densities are at 20 °C: 0.9982 g/cm^3 for H2O and 1.1053 g/cm^3 for D2O.
        constexpr double mass_D = 2.0141;
        constexpr double density_H2O = (2*scattering_length::H + scattering_length::O)*number_density(
            constants::mass::density::water, 2*constants::mass::get_mass(constants::atom_t::H) + constants::mass::get_mass(constants::atom_t::O)
        );
        constexpr double density_D2O = (2*scattering_length::D + scattering_length::O)*number_density(
            1.1053*constants::SI::mass::gm/constants::SI::volume::cm3, 2*mass_D + constants::mass::get_mass(constants::atom_t::O)
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
