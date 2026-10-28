// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/ParameterGenerationStrategy.h>

#include <math/Vector3.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <utility/Random.h>

#include <cassert>
#include <cmath>
#include <numbers>
#include <random>

using namespace ausaxs;
using namespace ausaxs::rigidbody::parameter;

ParameterGenerationStrategy::ParameterGenerationStrategy(
    observer_ptr<const Rigidbody> rigidbody, int iterations, const ParameterAmplitudes& amplitudes) 
    : rigidbody(rigidbody), amplitudes(amplitudes), decay_strategy(rigidbody::factory::create_decay_strategy(iterations, settings::rigidbody::DecayStrategyChoice::Linear)
) {}

ParameterGenerationStrategy::ParameterGenerationStrategy(
    observer_ptr<const Rigidbody> rigidbody, std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy, const ParameterAmplitudes& amplitudes) 
    : rigidbody(rigidbody), amplitudes(amplitudes), decay_strategy(std::move(decay_strategy)
) {}

ParameterGenerationStrategy::~ParameterGenerationStrategy() = default;

double ParameterGenerationStrategy::draw(double amplitude) {
    return unit_dist(random::generator())*amplitude;
}

Vector3<double> ParameterGenerationStrategy::draw_direction() {
    // inverse transform sampling
    double z = unit_dist(random::generator());
    double azimuth = std::numbers::pi*unit_dist(random::generator());
    double r = std::sqrt(1 - z*z);
    return {r*std::cos(azimuth), r*std::sin(azimuth), z};
}

Vector3<double> ParameterGenerationStrategy::draw_isotropic(double magnitude, int dimensions) {
    // a symmetry with no parameters of its own, e.g. a view onto another body's, has zero dimensions; a planar one has the first two
    assert(dimensions <= 3 && "ParameterGenerationStrategy::draw_isotropic: unexpected parameter count.");
    return draw_isotropic(magnitude, Axes{.x = 0 < dimensions, .y = 1 < dimensions, .z = 2 < dimensions});
}

Vector3<double> ParameterGenerationStrategy::draw_isotropic(double magnitude, const Axes& axes) {
    // a uniform direction projected onto a subspace is still uniform in direction within it, so only its length must be restored
    auto direction = draw_direction();
    if (!axes.x) {direction.x() = 0;}
    if (!axes.y) {direction.y() = 0;}
    if (!axes.z) {direction.z() = 0;}
    double length = direction.magnitude();
    return length == 0 ? Vector3<double>{0, 0, 0} : direction*(magnitude/length);
}

void ParameterGenerationStrategy::set_decay_strategy(std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy) {
    this->decay_strategy = std::move(decay_strategy);
}

observer_ptr<rigidbody::parameter::decay::DecayStrategy> ParameterGenerationStrategy::get_decay_strategy() const {return decay_strategy.get();}

void ParameterGenerationStrategy::set_max_translation_distance(double distance) {
    amplitudes.translation = distance;
}

void ParameterGenerationStrategy::set_max_rotation_angle(double radians) {
    amplitudes.rotation = radians;
}

void ParameterGenerationStrategy::set_max_symmetry_translation_distance(double distance) {
    amplitudes.symmetry_translation = distance;
}

void ParameterGenerationStrategy::set_max_symmetry_rotation_angle(double radians) {
    amplitudes.symmetry_rotation = radians;
}

const ParameterAmplitudes& ParameterGenerationStrategy::get_amplitudes() const {return amplitudes;}
