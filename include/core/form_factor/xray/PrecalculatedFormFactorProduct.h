// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ContainerFwd.h>
#include <constants/Constants.h>
#include <container/ArrayContainer2D.h>
#include <form_factor/xray/FormFactor.h>
#include <form_factor/FormFactorConcepts.h>

namespace ausaxs::form_factor {
    class PrecalculatedFormFactorProduct {
        public:
            constexpr PrecalculatedFormFactorProduct() noexcept = default;

            template<FormFactorType T1, FormFactorType T2>
            constexpr PrecalculatedFormFactorProduct(const T1& ff1, const T2& ff2) noexcept {
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

    namespace storage::atomic {
        using table_t = container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count(), form_factor::get_count()>;

        /**
         * @brief Get a precalculated atomic form factor product for a given pair of atomic form factors.
         * 
         * @param i The index of the first atomic form factor.
         * @param j The index of the second atomic form factor.
         */
        const PrecalculatedFormFactorProduct& get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept;

        /**
         * @brief Get the precalculated atomic form factor product table.
         *        The table is symmetric. 
         */
        const table_t& get_precalculated_form_factor_table() noexcept;
    }
}