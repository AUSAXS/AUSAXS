// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/neutron/NeutronFormFactorTable.h>
#include <math/ConstexprMath.h>
#include <utility/Exceptions.h>

#include <array>
#include <cmath>
#include <string>
#include <type_traits>
#include <utility>

namespace ausaxs::form_factor::neutron {
    /**
     * @brief The orientationally averaged neutron form factor of a heavy atom X with n equivalent bound hydrogens.
     *        Nuclei scatter as points, so both quantities below are exact for a rigid group and involve no parametrization.
     *        Scattering lengths are in fm, and may be negative.
     *
     *        The amplitude enters every pair of distinct groups, which contributes f_i(q)*f_j(q)*sinc(q*r_ij):
     *            f(q) = b_X + n*b_H*sinc(q*d_XH)
     *
     *        The self-term replaces f(q)^2 for a group with itself, since it retains the X-H and H-H interference within the group:
     *            s(q) = b_X^2 + n*b_H^2 + 2n*b_X*b_H*sinc(q*d_XH) + n(n-1)*b_H^2*sinc(q*d_HH)
     */
    class FormFactor {
        public:
            constexpr FormFactor() = default;

            /**
             * @param b_heavy The coherent scattering length of the heavy atom.
             * @param b_hydrogen The coherent scattering length of a single bound hydrogen.
             * @param hydrogens The number of bound hydrogens.
             * @param xh_distance The distance between the heavy atom and each hydrogen in Å.
             * @param hh_distance The distance between each pair of hydrogens in Å.
             */
            constexpr FormFactor(double b_heavy, double b_hydrogen, int hydrogens, double xh_distance, double hh_distance)
                : b_heavy(b_heavy), b_hydrogen(b_hydrogen), hydrogens(hydrogens), xh_distance(xh_distance), hh_distance(hh_distance) {}

            /**
             * @brief Initialize a form factor from a tabulated group.
             *
             * @param b_hydrogen The coherent scattering length of a single bound hydrogen.
             */
            constexpr FormFactor(const coefficients::Group& group, double b_hydrogen)
                : FormFactor(group.b_heavy, b_hydrogen, group.hydrogens, group.d_xh, group.d_hh) {}

            /**
             * @brief Evaluate the amplitude at a given q value.
             */
            constexpr double evaluate(double q) const {
                return b_heavy + hydrogens*b_hydrogen*sinc(q*xh_distance);
            }

            /**
             * @brief Evaluate the self-term at a given q value.
             */
            constexpr double evaluate_self(double q) const {
                return b_heavy*b_heavy + hydrogens*b_hydrogen*b_hydrogen
                    + 2*hydrogens*b_heavy*b_hydrogen*sinc(q*xh_distance)
                    + hydrogens*(hydrogens-1)*b_hydrogen*b_hydrogen*sinc(q*hh_distance);
            }

            /**
             * @brief Evaluate the amplitude at q = 0.
             */
            constexpr double I0() const {
                return b_heavy + hydrogens*b_hydrogen;
            }

        private:
            double b_heavy = 0;
            double b_hydrogen = 0;
            int hydrogens = 0;
            double xh_distance = 0;
            double hh_distance = 0;

            static constexpr double sinc(double x) {
                if (x == 0) {return 1;}
                if (std::is_constant_evaluated()) {return constexpr_math::sin(x)/x;}
                return std::sin(x)/x;
            }
    };

    namespace detail {
        // all hydrogens, implicit and explicit, are protium
        constexpr auto table = [] <std::size_t... I> (std::index_sequence<I...>) {
            return std::array<FormFactor, total_ff_count>{FormFactor(form_factor::detail::ff_info_table[I].neutron_coefficients, scattering_length::H)...};
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
}
