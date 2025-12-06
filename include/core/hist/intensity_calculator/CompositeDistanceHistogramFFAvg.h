// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief An instantiation of the CompositeDistanceHistogramFFAvgBase class that uses the default precalculated form factor table.
     *        For more information, see CompositeDistanceHistogramFFAvgBase.
     */
    class CompositeDistanceHistogramFFAvg : public CompositeDistanceHistogramFFAvgBase<form_factor::lookup::atomic::table_t> {
        using CompositeDistanceHistogramFFAvgBase::CompositeDistanceHistogramFFAvgBase;

        protected:
            const form_factor::lookup::atomic::table_t& get_ff_table() const override {
                return form_factor::lookup::atomic::raw::get_table();
            }
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFAvg>, "CompositeDistanceHistogramAvg should support nothrow move semantics.");
}