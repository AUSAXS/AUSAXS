#pragma once

#include <utility/Concepts.h>

#include <tuple>
#include <atomic>
#include <random>

template<numeric T> class Vector3;
namespace rigidbody {
    struct Parameter;
    
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
            ParameterGenerationStrategy(int iterations, double length_start, double rad_start);

            /**
             * @brief Destructor.
             */
            virtual ~ParameterGenerationStrategy();

            Parameter next();

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
}