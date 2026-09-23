// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsFwd.h>
#include <data/DataFwd.h>
#include <form_factor/ExvFormFactor.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/FormFactorType.h>
#include <math/ConstexprMath.h>

#include <array>
#include <cmath>
#include <utility>

namespace ausaxs::form_factor {
    class FormFactor {
        public:
            /**
             * @brief Initialize a vacuum form factor based on a 5-Gaussian approximation.
             */
            constexpr FormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : a(a), b(b), c(c) {}

            /**
             * @brief Initialize a vacuum form factor from a set of tabulated five-Gaussian coefficients.
             */
            constexpr FormFactor(const constants::form_factor::FiveGaussian& coefficients) : FormFactor(coefficients.a, coefficients.b, coefficients.c) {}

            /**
             * @brief Initialize an excluded volume form factor.
             *        This is only used to instantiate the average excluded volume form factor.
             *        Note that these excluded volume form factors are not normalized. 
             */
            constexpr FormFactor(const ExvFormFactor& ffx) : a({ffx.q0, 0, 0, 0, 0}), b({ffx.exponent, 0, 0, 0, 0}), c(0) {}

            /**
             * @brief Evaluate the form factor at a given q value.
             *        The vacuum form factors are normalized to 1 at q = 0.
             */
            constexpr double evaluate(double q) const {
                double sum = 0;
                for (int i = 0; i < 5; ++i) {
                    if (std::is_constant_evaluated()) {
                        sum += a[i]*constexpr_math::exp(-b[i]*q*q);
                    } else {
                        sum += a[i]*std::exp(-b[i]*q*q);
                    }
                }
                return (sum + c)*q0;
            }

            /**
             * @brief Evaluate the form factor at q = 0.
             */
            constexpr double I0() const {
                return (a[0] + a[1] + a[2] + a[3] + a[4] + c);
            }

            /**
             * @brief Manually set the normalization of this form factor.
             *        evaluate(0) will return this value.
             */
            constexpr void set_normalization(double q0) {
                this->q0 = q0/I0();
            }

        protected:
            double q0 = 1;

        private: 
            std::array<double, 5> a;
            std::array<double, 5> b;
            double c;
    };

    /**
     * The vacuum form factors of all form factor types, as described by form_factor::detail::ff_info_table.
     */
    namespace lookup::atomic::raw {
        namespace detail {
            constexpr auto table = [] <std::size_t... I> (std::index_sequence<I...>) {
                return std::array<FormFactor, total_ff_count>{FormFactor(form_factor::detail::ff_info_table[I].coefficients)...};
            }(std::make_index_sequence<total_ff_count>{});
        }

        constexpr const FormFactor& get(form_factor_t type) {
            if (type == form_factor_t::UNKNOWN) {
                throw ausaxs::except::runtime_error(
                    "form_factor::lookup::atomic::raw::get: Attempted to get the form factor of an UNKNOWN atom.\n"
                    "This typically occurs when performing species-dependent operations on data without form factor information."
                );
            }
            if (!form_factor::detail::is_tabulated(type)) {
                throw ausaxs::except::runtime_error("form_factor::lookup::atomic::raw::get: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
            }
            return detail::table[static_cast<int>(type)];
        }
    }
}