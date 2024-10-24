#pragma once

#include <vector>

namespace ausaxs::math {
    /**
     * @brief Find the indices of minima in the dataset.
     * 
     * @param x The x values of the dataset.
     * @param y The y values of the dataset.
     * @param min_spacing The minimum spacing between minima.
     * @param prominence The minimum prominence of a minima. This is the estimated depth of the minima. 
     */
    std::vector<unsigned int> find_minima(const std::vector<double>& x, const std::vector<double>& y, unsigned int min_spacing, double min_prominence);

    namespace detail {
        constexpr double min_slope = 1; // each point must be at least this much higher (in percent) than the previous point to be considered part of the bounds of a minima
    }
}