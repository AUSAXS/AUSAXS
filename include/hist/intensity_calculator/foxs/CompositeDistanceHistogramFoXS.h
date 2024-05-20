#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicitBase.h>
#include <hist/intensity_calculator/foxs/FormFactorFoXS.h>
#include <form_factor/PrecalculatedExvFormFactorProduct.h>

namespace hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that uses the same form factors as FoXS. 
     */
    class CompositeDistanceHistogramFoXS : public CompositeDistanceHistogramFFExplicitBase<form_factor::storage::atomic::table_t, form_factor::storage::cross::table_t, form_factor::storage::exv::table_t> {
        public: 
            using CompositeDistanceHistogramFFExplicitBase::CompositeDistanceHistogramFFExplicitBase;
            ~CompositeDistanceHistogramFoXS() override = default;

            const form_factor::storage::atomic::table_t& get_ff_table() const override {
                return ffaa_table;
            }

            const form_factor::storage::cross::table_t& get_ffax_table() const override {
                return ffax_table;
            }

            const form_factor::storage::exv::table_t& get_ffxx_table() const override {
                return ffxx_table;
            }

            Limit get_excluded_volume_scaling_factor_limits() const override;

        protected:
            double exv_factor(double q) const override;
            form_factor::storage::atomic::table_t ffaa_table = form_factor::foxs::storage::atomic::generate_table();
            form_factor::storage::cross::table_t ffax_table = form_factor::foxs::storage::cross::generate_table();
            form_factor::storage::exv::table_t ffxx_table = form_factor::foxs::storage::exv::generate_table();
    };
}