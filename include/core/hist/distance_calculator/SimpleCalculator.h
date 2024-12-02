#pragma once

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinates.h>

#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Simple interface to queue histogram calculations. 
     *        Submit the data and then call calculate to get the result.
     *        The caller must guarantee the lifetime of all submitted data.
     */
    template<bool weighted_bins>
    class SimpleCalculator {
        using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
        struct run_result {
            run_result(int size_self, int size_cross) : self(size_self), cross(size_cross) {}
            std::vector<GenericDistribution1D_t> self;
            std::vector<GenericDistribution1D_t> cross;
        };

        public:
            /**
             * @brief Queue a self-correlation calculation. 
             *        This is faster than calling the cross-correlation method with the same data, as some optimizations can be made. 
             *
             * @param a The data to calculate the self-correlation for. The reference must be valid until calculate is called.
             * @return The index of the data in the result vector.
             */
            int enqueue_calculate_self(const hist::detail::CompactCoordinates& a);

            /**
             * @brief Queue a cross-correlation calculation. 
             *
             * @param a1 The first set of data to calculate the cross-correlation for. The reference must be valid until calculate is called.
             * @param a2 The second set of data to calculate the cross-correlation for. The reference must be valid until calculate is called.
             * @return The index of the data in the result vector.
             */
            int enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2);

            /**
             * @brief Calculate the queued histograms. 
             *
             * @return The calculated histograms. 
             */
            run_result run();

        private:
            std::vector<std::reference_wrapper<const hist::detail::CompactCoordinates>> self;
            std::vector<std::reference_wrapper<const hist::detail::CompactCoordinates>> cross_1, cross_2;
    };
}