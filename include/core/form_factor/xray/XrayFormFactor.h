// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsFwd.h>
#include <data/DataFwd.h>
#include <form_factor/ExvFormFactor.h>
#include <form_factor/FormFactorInfo.h>
#include <form_factor/GroupFormFactor.h>
#include <form_factor/xray/XrayFormFactorTable.h>
#include <math/ConstexprMath.h>

#include <array>
#include <utility>

namespace ausaxs::form_factor::xray {
    /**
     * @brief The X-ray form factor of a single atom, described by a five-Gaussian approximation.
     */
    class AtomicFormFactor {
        public:
            /**
             * @brief Initialize a vacuum form factor based on a 5-Gaussian approximation.
             */
            constexpr AtomicFormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : a(a), b(b), c(c) {}

            /**
             * @brief Initialize a vacuum form factor from a set of tabulated five-Gaussian coefficients.
             */
            constexpr AtomicFormFactor(const coefficients::FiveGaussian& coefficients) : AtomicFormFactor(coefficients.a, coefficients.b, coefficients.c) {}

            /**
             * @brief Initialize an excluded volume form factor.
             *        This is only used to instantiate the average excluded volume form factor.
             *        Note that these excluded volume form factors are not normalized. 
             */
            constexpr AtomicFormFactor(const form_factor::detail::ExvFormFactor& ffx) : a({ffx.q0, 0, 0, 0, 0}), b({ffx.exponent, 0, 0, 0, 0}), c(0) {}

            /**
             * @brief Evaluate the form factor at a given q value.
             */
            constexpr double evaluate(double q) const {
                double sum = 0;
                for (int i = 0; i < 5; ++i) {
                    sum += a[i]*math::exp(-b[i]*q*q);
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
     * @brief The vacuum X-ray form factor of a form factor type: a heavy atom with its bound hydrogens.
     */
    using FormFactor = GroupFormFactor<AtomicFormFactor>;

    /**
     * The vacuum X-ray form factors of all form factor types, as described by form_factor::detail::ff_info_table.
     */
    namespace raw {
        namespace detail {
            // the free hydrogen atom, centred on its nucleus
            constexpr AtomicFormFactor hydrogen(coefficients::H);

            constexpr auto table = [] <std::size_t... I> (std::index_sequence<I...>) {
                return std::array<FormFactor, total_ff_count>{FormFactor(
                    AtomicFormFactor(form_factor::detail::ff_info_table[I].xray_coefficients), 
                    hydrogen, 
                    form_factor::detail::ff_info_table[I].geometry
                )...};
            }(std::make_index_sequence<total_ff_count>{});
        }

        constexpr const FormFactor& get(form_factor_t type) {
            if (type == form_factor_t::UNKNOWN) {
                throw ausaxs::except::runtime_error(
                    "form_factor::xray::raw::get: Attempted to get the form factor of an UNKNOWN atom.\n"
                    "This typically occurs when performing species-dependent operations on data without form factor information."
                );
            }
            if (!form_factor::detail::is_tabulated(type)) {
                throw ausaxs::except::runtime_error("form_factor::xray::raw::get: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
            }
            return detail::table[static_cast<int>(type)];
        }
    }
}