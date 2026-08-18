// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/Rigidbody.h>
#include <math/Vector3.h>
#include <utility/Random.h>

#include <numbers>
#include <random>

using namespace ausaxs;
using namespace ausaxs::rigidbody::parameter;

double rigidbody::parameter::default_rotation() {
    return std::numbers::pi/3;
}

double rigidbody::parameter::default_symmetry_translation(observer_ptr<const Rigidbody> rigidbody) {
    return 2*rigidbody->molecule.get_Rg(false);
}

ParameterAmplitudes rigidbody::parameter::default_amplitudes(observer_ptr<const Rigidbody> rigidbody) {
    return {
        .translation = default_translation(),
        .rotation = default_rotation(),
        .symmetry_translation = default_symmetry_translation(rigidbody),
        .symmetry_rotation = default_symmetry_rotation()
    };
}

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
