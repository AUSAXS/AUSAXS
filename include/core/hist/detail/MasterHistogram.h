// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/GenericDistribution1D.h>
#include <hist/Histogram.h>
#include <utility/Axis.h>

namespace ausaxs::hist::detail {
    // Simple typedef for clarity.
    template<bool use_weighted_distribution>
    using PartialHistogram = typename GenericDistribution1D<use_weighted_distribution>::type;

    // Simple typedef for clarity.
    template<bool use_weighted_distribution>
    using HydrationHistogram = typename GenericDistribution1D<use_weighted_distribution>::type;

    /**
     * @brief We also define the MasterHistogram type, which is identical to a PartialHistogram. 
     *        We do this to make += and -= well-defined operations. 
     */
    template<bool use_weighted_distribution>
    class MasterHistogram : public GenericDistribution1D<use_weighted_distribution>::type {
        using GenericDistribution1D_t = typename GenericDistribution1D<use_weighted_distribution>::type;
        public: 
            MasterHistogram();

            /**
             * @brief Create a new Master Histogram. 
             * @param p The current histogram. 
             * @param p_base The constant, unchanging part of the histogram. 
             */
            MasterHistogram(const std::vector<double>& p_base, const Axis& axis);

            /**
             * @brief Create a new Master Histogram. 
             * @param p The current histogram. 
             * @param p_base The constant, unchanging part of the histogram. 
             */
            MasterHistogram(std::vector<double>&& p_base, const Axis& axis);

            /**
             * @brief Add a PartialHistogram to the MasterHistogram. 
             */
            MasterHistogram& operator+=(const GenericDistribution1D_t& rhs);

            /**
             * @brief Subtract a PartialHistogram from the MasterHistogram. We have to use a lambda since the standard std::minus would
             *        reverse the order of the entries.
             */
            MasterHistogram& operator-=(const GenericDistribution1D_t& rhs);

            // The base part of the histogram which will never change. This contains all internal distances between atoms in each individual body.
            GenericDistribution1D_t base; //? remove?
            Axis axis;
    };
}