#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief An instantiation of the CompositeDistanceHistogramFFAvgBase class that uses the default precalculated form factor table.
     *        For more information, see CompositeDistanceHistogramFFAvgBase.
     */
    class CompositeDistanceHistogramFFAvg : public CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t> {
        using CompositeDistanceHistogramFFAvgBase::CompositeDistanceHistogramFFAvgBase;

        protected:
            const form_factor::storage::atomic::table_t& get_ff_table() const override {
                return form_factor::storage::atomic::get_precalculated_form_factor_table();
            }
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFAvg>, "CompositeDistanceHistogramAvg should support nothrow move semantics.");
}