// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/FormFactorLookupFwd.h>
#include <form_factor/FormFactorType.h>
#include <table/DebyeTableManager.h>
#include <utility/observer_ptr.h>

#include <vector>

namespace ausaxs::hist {
    /**
     * @brief Common base for the excluded volume models built on a space-filling grid of dummy spheres.
     */
    class CompositeDistanceHistogramFFGridBase : public CompositeDistanceHistogramFFAvgBase<form_factor::lookup::table_t> {
        public:
            using CompositeDistanceHistogramFFAvgBase::CompositeDistanceHistogramFFAvgBase;

            const form_factor::lookup::table_t& get_ff_table() const override {return ff_table;}

            /**
             * @brief Regenerate the internal form factor table for the grid-based calculations.
             *        This must be called to reflect changes in settings::grid::exv::width.
             * @param ffx The excluded volume form factor to use.
             */
            template<FormFactorType T>
            static void regenerate_ff_table(T&& ffx);

            /**
             * @brief Regenerate the internal form factor table, coupling the excluded volume form factor to the grid volume.
             */
            static void regenerate_ff_table();

            /**
             * @brief Get the distance axis for the excluded volume calculations.
             *        If weighted bins are used, this will be distinct from the regular distance axis.
             */
            const std::vector<double>& get_d_axis_xx() const {return distance_axes.xx;}

            /**
             * @brief Get the distance axis for the cross term calculations.
             *        If weighted bins are used, this will be distinct from the regular distance axis.
             */
            const std::vector<double>& get_d_axis_ax() const {return distance_axes.ax;}

        protected:
            /**
             * @brief Get the sinc(x) lookup table for the excluded volume for the Debye transform.
             */
            observer_ptr<const table::DebyeTable> get_sinc_table_xx() const;

            /**
             * @brief Get the sinc(x) lookup table for the cross terms for the Debye transform.
             */
            observer_ptr<const table::DebyeTable> get_sinc_table_ax() const;

            /**
             * @brief Store the grid distance axes and build their sinc(qd) lookup tables.
             */
            void initialize_grid_axes(std::vector<double>&& d_axis_ax, std::vector<double>&& d_axis_xx);

            inline static form_factor::lookup::table_t ff_table;
            struct {table::DebyeTableManager xx, ax;} sinc_tables;
            struct {std::vector<double> xx, ax;} distance_axes;
    };
}
