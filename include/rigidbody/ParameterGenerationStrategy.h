#pragma once

#include <tuple>
#include <random>

#include "math/Vector3.h"
#include "rigidbody/Parameters.h"

class ParameterGenerationStrategy {
  public:
    /**
     * @brief Constructor.
     * 
     * @param iterations The expected number of iterations. 
     * @param length_start The start length of the generated translation vectors. 
     * @param rad_start The start angle in radians of the generated rotations. 
     */
    ParameterGenerationStrategy(const int& iterations, const double& length_start, const double& rad_start) : iterations(iterations) {
      std::random_device random;
      generator = std::mt19937(random());
      translation_dist = std::uniform_int_distribution<int>(-length_start, length_start);
      rotation_dist = std::uniform_int_distribution<int>(-rad_start, rad_start);
    }

    /**
     * @brief Destructor.
     */
    virtual ~ParameterGenerationStrategy() = default;

    Parameters::Parameter next() {
      auto[rx, ry, rz] = get_rotation();
      Vector3 x = get_translation();
      iteration++;
      return Parameters::Parameter(x, rx, ry, rz);
    }

  protected:
    int iteration = 0;                                   // Current iteration. 
    int iterations;                                      // The total number of iterations. Used to determine the current scaling. 
    std::mt19937 generator;                              // The random number generator. 
    std::uniform_int_distribution<int> translation_dist; // The random number distribution for translations. 
    std::uniform_int_distribution<int> rotation_dist;    // The random number distribution for rotations. 

    /**
     * @brief Get a new rotation offset based on the current iteration. 
     */
    virtual std::tuple<double, double, double> get_rotation() = 0;

    /**
     * @brief Get a new translation offset based on the current iteration. 
     */
    virtual Vector3 get_translation() = 0;
};