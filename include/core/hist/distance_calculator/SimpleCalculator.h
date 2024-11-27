#pragma once

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/detail/CompactCoordinates.h>

#include <vector>

namespace ausaxs::hist::detail {
    /**
     * @brief Simple interface to queue histogram calculations. 
     *        Submit the data and then call calculate to get the result.
     *        The caller must guarantee the lifetime of all submitted data.
     */
    class SimpleCalculator {
        public:
            void enqueue_calculate_self(const hist::detail::CompactCoordinates& a);
            void enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2);
            std::unique_ptr<ICompositeDistanceHistogram> run();

        private:
            std::vector<std::reference_wrapper<const hist::detail::CompactCoordinates>> self;
            std::vector<std::reference_wrapper<const hist::detail::CompactCoordinates>> cross_1, cross_2;
    };
}