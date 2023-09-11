#pragma once

#include <hist/Histogram.h>

#include <vector>

namespace hist::detail {
    /**
     * @brief Simple data containers defined for clarity.  
     */
    typedef Histogram PartialHistogram;
    typedef Histogram HydrationHistogram;

    /**
     * @brief We also define the MasterHistogram type, which is identical to a PartialHistogram. 
     *        We do this to make += and -= well-defined operations. 
     */
    class MasterHistogram : public Histogram {
        public: 
            MasterHistogram() {}

            /**
             * @brief Create a new Master Histogram. 
             * @param p The current histogram. 
             * @param p_base The constant, unchanging part of the histogram. 
             */
            MasterHistogram(const std::vector<double>& p_base, const Axis& axis);

            /**
             * @brief Add a PartialHistogram to the MasterHistogram. 
             */
            MasterHistogram& operator+=(const PartialHistogram& rhs);

            /**
             * @brief Subtract a PartialHistogram from the MasterHistogram. We have to use a lambda since the standard std::minus would
             *        reverse the order of the entries.
             */
            MasterHistogram& operator-=(const PartialHistogram& rhs);

            // The base part of the histogram which will never change. This contains all internal distances between atoms in each individual body.
            Histogram base;
    };
}