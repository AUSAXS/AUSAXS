/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/parameters/Parameters.h>
#include <math/Vector3.h>

#include <random>

using namespace rigidbody::parameter;

ParameterGenerationStrategy::ParameterGenerationStrategy(int iterations, double length_start, double rad_start) : iterations(iterations) {
    std::random_device random;
    generator = std::mt19937(random());
    translation_dist = std::uniform_real_distribution<double>(-length_start, length_start);
    rotation_dist = std::uniform_real_distribution<double>(-rad_start, rad_start);
}

ParameterGenerationStrategy::~ParameterGenerationStrategy() = default;

Parameter ParameterGenerationStrategy::next() {
    auto[rx, ry, rz] = get_rotation();
    Vector3 x = get_translation();
    iteration++;
    return Parameter(x, rx, ry, rz);
}
