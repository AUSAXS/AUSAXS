#include <rigidbody/parameters/ParameterGenerationStrategy.h>

using namespace rigidbody;

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