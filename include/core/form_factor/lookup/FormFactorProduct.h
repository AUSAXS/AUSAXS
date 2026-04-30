// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/FormFactorConcepts.h>
#include <container/ArrayContainer2D.h>

namespace ausaxs::form_factor {
    class FormFactorProduct {
        public:
            constexpr FormFactorProduct() noexcept = default;

            template<FormFactorType T1, FormFactorType T2>
            constexpr FormFactorProduct(const T1& ff1, const T2& ff2) noexcept {
                std::array<double, constants::axes::q_axis.bins> res;
                for (unsigned int i = 0; i < res.size(); ++i) {
                    res[i] = ff1.evaluate(constants::axes::q_vals[i])*ff2.evaluate(constants::axes::q_vals[i]);
                }
                precalculated_ff_q = std::move(res);
            }

            /**
             * @brief Get the precalculated form factor product for a given q value.
             * 
             * These products are calculated at compile-time for the default q axis defined in the constants namespace.
             */
            constexpr double evaluate(unsigned int index) const noexcept {
                return precalculated_ff_q[index];
            }

        private:
            std::array<double, constants::axes::q_axis.bins> precalculated_ff_q;
    };
}