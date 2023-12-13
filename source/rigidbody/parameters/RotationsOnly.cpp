#include <rigidbody/parameters/RotationsOnly.h>
#include <math/Vector3.h>

#include <random>

using namespace rigidbody::parameter;

RotationsOnly::RotationsOnly(int iterations, double length_start, double rad_start) : ParameterGenerationStrategy(iterations, length_start, rad_start) {}

RotationsOnly::~RotationsOnly() = default;

std::tuple<double, double, double> RotationsOnly::get_rotation() {
    double scaling = scale();

    double dr1 = rotation_dist(generator)*scaling;
    double dr2 = rotation_dist(generator)*scaling;
    double dr3 = rotation_dist(generator)*scaling;
    return std::tuple(dr1, dr2, dr3);
}

Vector3<double> RotationsOnly::get_translation() {
    return Vector3(0., 0., 0.);
}

double RotationsOnly::scale() const {
    return double(iterations - iteration)/iterations; 
}