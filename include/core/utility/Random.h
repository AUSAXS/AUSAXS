// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/Logging.h>

#include <random>

namespace ausaxs::random {
    namespace detail {
        inline std::mt19937& get_generator() {
            static std::mt19937 gen = [](){
                int seed = std::random_device{}();
                logging::log("RNG: Initializing random number generator with seed: " + std::to_string(seed));
                return std::mt19937(seed);
            }();
            return gen;
        } 
    }

    /**
     * @brief Set the seed for the random number generator.
     *        A new Mersenne Twister generator with this seed will be initialized and used throughout the library.
     */
    inline void set_seed(int seed) noexcept {
        detail::get_generator().seed(seed);
        logging::log("RNG: Setting random number generator seed to: " + std::to_string(seed));
    }

    inline std::mt19937& generator() {
        return detail::get_generator();
    }
}