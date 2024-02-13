#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>

namespace hist {
    /**
     * @brief An instantiation of the CompositeDistanceHistogramFFAvgBase class that uses the default precalculated form factor table.
     *        For more information, see CompositeDistanceHistogramFFAvgBase.
     */
    class CompositeDistanceHistogramFFAvg : public CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t> {
        using CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>::CompositeDistanceHistogramFFAvgBase;

        const form_factor::storage::atomic::table_t& get_ff_table() const override {
            return form_factor::storage::atomic::get_precalculated_form_factor_table();
        }
    };
}