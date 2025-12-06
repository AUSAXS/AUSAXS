// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicitBase.h>
#include <form_factor/lookup/ExvFormFactorProduct.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions and atomic types. 
     *        Beyond the functionality of CompositeDistanceHistogram, this class also uses individual form factors for each atomic type.
     *        Unique exluded volume form factors are used for each atomic type.
     *        For more information, see CompositeDistanceHistogram.
     */
    class CompositeDistanceHistogramFFExplicit : public CompositeDistanceHistogramFFExplicitBase<form_factor::lookup::atomic::table_t, form_factor::lookup::cross::table_t, form_factor::lookup::exv::table_t> {
        public: 
            using CompositeDistanceHistogramFFExplicitBase::CompositeDistanceHistogramFFExplicitBase;

            const form_factor::lookup::atomic::table_t& get_ff_table() const override {
                return form_factor::lookup::atomic::raw::get_table();
            }

            const form_factor::lookup::cross::table_t& get_ffax_table() const override {
                return form_factor::lookup::cross::raw::get_table();
            }

            const form_factor::lookup::exv::table_t& get_ffxx_table() const override {
                return form_factor::lookup::exv::raw::get_table();
            }
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFExplicit>, "CompositeDistanceHistogramFFExplicit should be nothrow move constructible");
}