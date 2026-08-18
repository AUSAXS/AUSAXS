// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/Rigidbody.h>
#include <math/Vector3.h>
#include <utility/Random.h>

#include <cassert>
#include <cmath>
#include <numbers>
#include <random>

using namespace ausaxs;
using namespace ausaxs::rigidbody::parameter;

ParameterGenerationStrategy::ParameterGenerationStrategy(
    observer_ptr<const Rigidbody> rigidbody, unsigned int iterations, const ParameterAmplitudes& amplitudes) 
    : rigidbody(rigidbody), amplitudes(amplitudes), decay_strategy(rigidbody::factory::create_decay_strategy(iterations)
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

Vector3<double> ParameterGenerationStrategy::draw_isotropic(double magnitude, std::size_t dimensions) {
    switch (dimensions) {
        case 0: return {0, 0, 0}; // a symmetry with no parameters of its own, e.g. a view onto another body's
        case 2: { // a planar parameter, whose direction is a single uniform azimuth
            double azimuth = std::numbers::pi*unit_dist(random::generator());
            return {magnitude*std::cos(azimuth), magnitude*std::sin(azimuth), 0};
        }
        case 3: return draw_direction()*magnitude;
        default:
            assert(false && "ParameterGenerationStrategy::draw_isotropic: unexpected parameter count.");
            return {0, 0, 0};
    }
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
