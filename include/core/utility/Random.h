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

        inline std::normal_distribution<double>& get_gaussian() {
            static std::normal_distribution<double> dist;
            return dist;
        }
    }

    /**
     * @brief Set the seed for the random number generator.
     *        A new Mersenne Twister generator with this seed will be initialized and used throughout the library.
     */
    inline void set_seed(int seed) noexcept {
        detail::get_generator().seed(seed);

        // the gaussian generates its variates in pairs and keeps the unused one around between calls. that leftover would
        // survive the reseed and make the first draw depend on how many draws preceded it, so it has to be discarded.
        detail::get_gaussian().reset();
        logging::log("RNG: Setting random number generator seed to: " + std::to_string(seed));
    }

    inline std::mt19937& generator() {
        return detail::get_generator();
    }

    /**
     * @brief Draw a normally distributed number from the library-wide generator.
     *        Use this instead of a local std::normal_distribution, since only the shared one is reset by set_seed.
     */
    inline double gaussian(double mean, double stddev) {
        return detail::get_gaussian()(detail::get_generator(), std::normal_distribution<double>::param_type(mean, stddev));
    }
}