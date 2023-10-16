#pragma once

#include <container/ContainerFwd.h>
#include <constants/Constants.h>
#include <form_factor/FormFactor.h>
#include <container/ArrayContainer2D.h>

#include <vector>

namespace form_factor {
    class FormFactor;
    class PrecalculatedFormFactorProduct {
        public:
            constexpr PrecalculatedFormFactorProduct() noexcept = default;
            constexpr PrecalculatedFormFactorProduct(const FormFactor& ff1, const FormFactor& ff2) noexcept {
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

    namespace storage {
        const PrecalculatedFormFactorProduct& get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept;
        const container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count(), form_factor::get_count()>& get_precalculated_form_factor_table() noexcept;
    }
}