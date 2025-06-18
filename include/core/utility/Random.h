#pragma once

#include <random>

namespace ausaxs::random {
    namespace detail {
        inline std::mt19937& get_generator() {
            static std::mt19937 gen(std::random_device{}());
            return gen;
        } 
    }

    /**
     * @brief Set the seed for the random number generator.
     *        A new Mersenne Twister generator with this seed will be initialized and used throughout the library.
     */
    inline void set_seed(int seed) noexcept {
        detail::get_generator().seed(seed);
    }

    inline std::mt19937& generator() {
        return detail::get_generator();
    }
}