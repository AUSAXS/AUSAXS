#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/crysol/ExvFormFactorCrysol.h>

namespace hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that mimics the CRYSOL excluded volume fitting. 
     */
    
    class CompositeDistanceHistogramCrysol : public CompositeDistanceHistogramFFExplicitBase<form_factor::storage::atomic::table_t, form_factor::storage::cross::table_t, form_factor::storage::exv::table_t>{
        public:
            using CompositeDistanceHistogramFFExplicitBase::CompositeDistanceHistogramFFExplicitBase;
            ~CompositeDistanceHistogramCrysol() override = default;

            const form_factor::storage::atomic::table_t& get_ff_table() const override {
                return form_factor::storage::atomic::get_precalculated_form_factor_table();
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

            form_factor::storage::cross::table_t ffax_table = form_factor::crysol::storage::cross::generate_table();
            form_factor::storage::exv::table_t ffxx_table = form_factor::crysol::storage::exv::generate_table();
    };
}