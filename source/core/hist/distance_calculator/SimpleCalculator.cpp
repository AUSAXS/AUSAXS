#include <hist/distance_calculator/SimpleCalculator.h>

using namespace ausaxs;

void hist::detail::SimpleCalculator::enqueue_calculate_self(const hist::detail::CompactCoordinates& a) {
    self.push_back(a);
}

void hist::detail::SimpleCalculator::enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2) {
    cross_1.push_back(a1);
    cross_2.push_back(a2);
}

std::unique_ptr<hist::ICompositeDistanceHistogram> hist::detail::SimpleCalculator::run() {}