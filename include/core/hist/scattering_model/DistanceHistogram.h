// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>
#include <dataset/DatasetFwd.h>
#include <hist/HistFwd.h>
#include <hist/Histogram.h>
#include <hist/ScatteringProfile.h>
#include <hist/distribution/DistributionFwd.h>
#include <table/DebyeTableManager.h>
#include <utility/Axis.h>
#include <utility/TypeTraits.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <span>
#include <vector>

namespace ausaxs::hist {
    /**
     * @brief A class containing a single total distance histogram for use in the Debye transform.
     *        This is the simplest implementation of the transform; for more advanced features, see its subclasses.
     */
    class DistanceHistogram : protected Histogram {
        public: 
            DistanceHistogram();
            DistanceHistogram(const DistanceHistogram&) = delete;
            DistanceHistogram(DistanceHistogram&& other) noexcept;
            DistanceHistogram& operator=(const DistanceHistogram&) = delete;
            DistanceHistogram& operator=(DistanceHistogram&& other) noexcept;

            /**
             * @brief Create an unweighted distance histogram.
             */
            DistanceHistogram(hist::Distribution1D&& p_tot);

            /**
             * @brief Create a weighted distance histogram. 
             *        A custom sinc(x) lookup table based on the weights in @a p_tot will be calculated and used in the Debye transform.
             */
            DistanceHistogram(hist::WeightedDistribution1D&& p_tot);

            /**
             * @brief Create a distance histogram from a scattering model.
             *        Only the total distance histogram will be extracted from the scattering model.
             */
            DistanceHistogram(std::unique_ptr<IScatteringModel> model);

            ~DistanceHistogram() override;

            /**
             * @brief Perform the Fourier transform through the Debye equation.
             */
            virtual ScatteringProfile debye_transform() const;

            /**
             * @brief Perform the Fourier transform through the Debye equation.
             *
             * @param q The q values at which to evaluate the scattering.
             */
            Dataset debye_transform(const std::vector<double>& q) const;

            /**
             * @brief Get the distance axis describing the current histogram.
             */
            const std::vector<double>& get_d_axis() const;

            /**
             * @brief Get the q axis used in the Fourier transform. 
             */
            static const std::vector<double>& get_q_axis();

            /**
             * @brief Get the total histogram counts. Equivalent to get_counts().
             */
            virtual const std::vector<double>& get_weighted_counts() const;

            /**
             * @brief Determine if the structure represented by this histogram is highly ordered.
             *        This is done by estimating the number of peaks in the histogram and comparing it to the number of bins.             
             */
            bool is_highly_ordered() const;

            /**
             * @brief Determine whether a distance-count sequence is highly ordered.
             *        This is done by estimating the number of peaks in the histogram and comparing it to the number of bins.             
             */
            static bool is_highly_ordered(const std::vector<double>& counts);

        protected:
            /**
             * @brief The Debye transform I(q) = sum_d counts(d)*sinc(q*d) of a single count array, over the current q range.
             *
             * @param sinqd_table A sinc table over the full default q axis.
             */
            static ScatteringProfile transform_counts(std::span<const double> counts, observer_ptr<const table::DebyeTable> sinqd_table);

            /**
             * @brief As transform_counts, but at the given q values.
             *
             * @param sinqd_table The sinc table. Its row @a first_row + i must hold the sinc values for q[i].
             */
            static std::vector<double> transform_counts_at(
                std::span<const double> counts, observer_ptr<const table::DebyeTable> sinqd_table, std::span<const double> q, int first_row = 0
            );

            /**
             * @brief Evaluate the intensity directly at arbitrary q values, without using precomputed tables.
             */
            virtual std::vector<double> debye_transform_direct(const std::vector<double>& q) const;

            /**
             * @brief A sinc table for the distance axis of this histogram at the given q values, with one row per q value.
             */
            observer_ptr<const table::DebyeTable> sinc_table_for(const std::vector<double>& q) const;
            //! Unclear, change name/function. 

            /**
             * @brief Expand settings::axes::qmin..qmax to cover @a q, which must lie within the default q axis.
             *        Subclasses caching anything sized by the q range must invalidate it.
             */
            virtual void expand_q_range(const std::vector<double>& q) const;
            //! Should probably be a pair of `expand_q_range` and `virtual on_q_range_changed` instead. 

            std::vector<double> d_axis;                     // The distance axis.
            mutable table::DebyeTableManager sinc_table;    // The sinc(x) lookup table manager for the Debye transform.

        private:
            void initialize();
            void initialize(std::vector<double>&& d_axis);
    };
    static_assert(supports_nothrow_move_v<DistanceHistogram>, "DistanceHistogram should support nothrow move semantics.");
}