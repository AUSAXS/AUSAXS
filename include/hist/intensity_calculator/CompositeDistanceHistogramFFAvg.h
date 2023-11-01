#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>

namespace hist {
    template<bool use_weighted_distribution>
    class CompositeDistanceHistogramFFAvg : public CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t, use_weighted_distribution> {
        using CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t, use_weighted_distribution>::CompositeDistanceHistogramFFAvgBase;

        const form_factor::storage::atomic::table_t& get_ff_table() const override {
            return form_factor::storage::atomic::get_precalculated_form_factor_table();
        }
    };
}