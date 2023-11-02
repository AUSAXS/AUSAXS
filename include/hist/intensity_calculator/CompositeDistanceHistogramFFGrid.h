#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>

namespace hist {
    class CompositeDistanceHistogramFFGrid : public CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t> {
        using CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>::CompositeDistanceHistogramFFAvgBase;

        const form_factor::storage::atomic::table_t& get_ff_table() const override {
            static auto ff_table = generate_table();
            return ff_table;
        }

        static form_factor::storage::atomic::table_t generate_table();
    };
}