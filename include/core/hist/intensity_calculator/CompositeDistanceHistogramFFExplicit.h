#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicitBase.h>
#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions and atomic types. 
     *        Beyond the functionality of CompositeDistanceHistogram, this class also uses individual form factors for each atomic type.
     *        Unique exluded volume form factors are used for each atomic type.
     *        For more information, see CompositeDistanceHistogram.
     */
    class CompositeDistanceHistogramFFExplicit : public CompositeDistanceHistogramFFExplicitBase<form_factor::storage::atomic::table_t, form_factor::storage::cross::table_t, form_factor::storage::exv::table_t> {
        public: 
            using CompositeDistanceHistogramFFExplicitBase::CompositeDistanceHistogramFFExplicitBase;

            const form_factor::storage::atomic::table_t& get_ff_table() const override {
                return form_factor::storage::atomic::get_precalculated_form_factor_table();
            }

            const form_factor::storage::cross::table_t& get_ffax_table() const override {
                return form_factor::storage::cross::get_precalculated_form_factor_table();
            }

            const form_factor::storage::exv::table_t& get_ffxx_table() const override {
                return form_factor::storage::exv::get_precalculated_form_factor_table();
            }
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFExplicit>, "CompositeDistanceHistogramFFExplicit should be nothrow move constructible");
}