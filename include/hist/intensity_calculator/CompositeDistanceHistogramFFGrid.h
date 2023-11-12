#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>

namespace hist {
    /**
     * @brief A composite distance histogram that uses a grid-based approach to calculate the form factor.
     *        All excluded volume form factors are considered equivalent. 
     */
    class CompositeDistanceHistogramFFGrid : public CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t> {
        public:
            using CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>::CompositeDistanceHistogramFFAvgBase;

            /**
             * @brief Get the form factor table for the grid-based calculations.
             */
            const form_factor::storage::atomic::table_t& get_ff_table() const override {return ff_table;}

            /**
             * @brief Regenerate the form factor table.
             *        This is only necessary if the excluded volume radius has changed.
             */
            static void regenerate_table() {ff_table = generate_table();}

        private: 
            static form_factor::storage::atomic::table_t generate_table();
            inline static form_factor::storage::atomic::table_t ff_table = generate_table();
    };
}