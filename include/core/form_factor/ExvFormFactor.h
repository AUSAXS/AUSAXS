// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>
#include <form_factor/ExvTable.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/FormFactorType.h>
#include <math/ConstexprMath.h>

#include <array>
#include <cmath>
#include <numbers>
#include <optional>

namespace ausaxs::form_factor {
    /**
     * @brief Calculate the excluded volume form factor based on the description from Fraser, MacRae & Suzuki: https://doi.org/10.1107/S0021889878014296
     */
    class ExvFormFactor {
        public: 
            /**
             * @brief Create a new excluded volume form factor with the given volume.
             *
             * @param volume The excluded volume of the atom. 
             */
            constexpr ExvFormFactor(double volume) 
                : exponent(constexpr_math::pow(volume, 2./3)/(4*std::numbers::pi)), q0(volume*constants::charge::density::water) 
            {}

            constexpr double evaluate_normalized(double q) const {
                if (std::is_constant_evaluated()) {
                    return constexpr_math::exp(-exponent*q*q);
                }
                return std::exp(-exponent*q*q);
            }

            constexpr double evaluate(double q) const {
                return q0*evaluate_normalized(q);
            }

            constexpr bool is_initialized() const {
                return exponent != 0;
            }

            double exponent = 0;
            double q0 = 1;
    };

    namespace detail {
        /**
         * @brief The excluded volume form factors of a single displaced volume set. 
         *        Form factor types absent from the volume set are also absent here. 
         */
        struct ExvFormFactorSet {
            constexpr ExvFormFactorSet(const constants::exv::detail::ExvSet& set) {
                for (int i = 0; i < total_ff_count; ++i) {
                    if (set.volumes[i].has_value()) {form_factors[i] = ExvFormFactor(*set.volumes[i]);}
                }
            }

            /**
             * @brief Check if this set has an excluded volume form factor for the given form factor type.
             */
            constexpr bool contains(form_factor_t type) const {
                return form_factor::detail::is_tabulated(type) && form_factors[static_cast<int>(type)].has_value();
            }

            constexpr ExvFormFactor get(form_factor_t type) const {
                if (!contains(type)) {
                    throw ausaxs::except::runtime_error("form_factor::detail::ExvFormFactorSet::get: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
                }
                return *form_factors[static_cast<int>(type)];
            }

            std::array<std::optional<ExvFormFactor>, total_ff_count> form_factors;
        };
    }
}