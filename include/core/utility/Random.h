#pragma once

#include <random>

namespace ausaxs::random {
    /**
     * @brief Set the seed for the random number generator.
     *        A new Mersenne Twister generator with this seed will be initialized and used throughout the library.
     */
    void set_seed(int seed) noexcept;
    const std::mt19937& generator();
}