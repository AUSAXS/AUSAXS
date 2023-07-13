#pragma once

#include <vector>

namespace math {
    /**
     * @brief Find the indices of minima in the dataset.
     * 
     * @param x The x values of the dataset.
     * @param y The y values of the dataset.
     * @param min_spacing The minimum spacing between minima.
     * @param prominence The minimum prominence of a minima. This is the estimated depth of the minima. 
     */
    std::vector<unsigned int> find_minima(const std::vector<double>& x, std::vector<double> y, unsigned int min_spacing, double min_prominence);
}