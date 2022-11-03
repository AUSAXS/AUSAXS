#pragma once

#include <tuple>
#include <random>
#include <atomic>

#include <math/Vector3.h>
#include <rigidbody/Parameters.h>

/**
 * @brief Thread-safe superclass for parameter generation strategies.
 */
class ParameterGenerationStrategy {
    public:
        /**
         * @brief Constructor.
         * 
         * @param iterations The expected number of iterations. 
         * @param length_start The start length of the generated translation vectors. 
         * @param rad_start The start angle in radians of the generated rotations. 
         */
        ParameterGenerationStrategy(int iterations, double length_start, double rad_start) : iterations(iterations) {
          std::random_device random;
          generator = std::mt19937(random());
          translation_dist = std::uniform_real_distribution<double>(-length_start, length_start);
          rotation_dist = std::uniform_real_distribution<double>(-rad_start, rad_start);
        }

        /**
         * @brief Destructor.
         */
        virtual ~ParameterGenerationStrategy() = default;

        Parameter next() {
          auto[rx, ry, rz] = get_rotation();
          Vector3 x = get_translation();
          iteration++;
          return Parameter(x, rx, ry, rz);
        }

    protected:
        std::atomic_uint iteration = 0;                          // Current iteration. 
        int iterations;                                          // The total number of iterations. Used to determine the current scaling. 
        std::mt19937 generator;                                  // The random number generator. 
        std::uniform_real_distribution<double> translation_dist; // The random number distribution for translations. 
        std::uniform_real_distribution<double> rotation_dist;    // The random number distribution for rotations. 

        /**
         * @brief Get a new rotation offset based on the current iteration. 
         */
        virtual std::tuple<double, double, double> get_rotation() = 0;

        /**
         * @brief Get a new translation offset based on the current iteration. 
         */
        virtual Vector3<double> get_translation() = 0;
};