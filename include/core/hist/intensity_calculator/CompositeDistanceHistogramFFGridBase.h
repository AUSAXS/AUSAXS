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
     *
     * Unlike the averaged model (CompositeDistanceHistogramFFAvg), the grid excluded volume is computed
     * independently of the atomic positions, so it cannot be derived from the atomic distance data. It is
     * stored in the excluded volume row and column of the distance histograms, and - since the grid is
     * highly ordered - it needs its own weighted distance axes and sinc(qd) lookup tables.
     *
     * This class owns that machinery, along with the grid form factor table shared by all grid variants.
     */
    class CompositeDistanceHistogramFFGridBase : public CompositeDistanceHistogramFFAvgBase<form_factor::lookup::table_t> {
        public:
            using CompositeDistanceHistogramFFAvgBase::CompositeDistanceHistogramFFAvgBase;

            const form_factor::lookup::table_t& get_ff_table() const override {return ff_table;}

            /**
             * @brief Generate a new internal form factor table for the grid-based calculations.
             *
             * @param ffx The excluded volume form factor to use. Leave as default to couple it to the grid volume.
             */
            template<FormFactorType T>
            static form_factor::lookup::table_t generate_ff_table(T&& ffx);
            static form_factor::lookup::table_t generate_ff_table(); //< @copydoc generate_ff_table(T&&)

            /**
             * @brief Regenerate the form factor table. This must be called to reflect changes in settings::grid::exv::width.
             *
             * @param ffx The excluded volume form factor to use. Leave as default to couple it to the grid volume.
             */
            template<FormFactorType T>
            static void regenerate_ff_table(T&& ffx = {0});

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
