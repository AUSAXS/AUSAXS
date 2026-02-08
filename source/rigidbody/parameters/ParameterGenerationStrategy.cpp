// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/Rigidbody.h>
#include <math/Vector3.h>
#include <utility/Random.h>

#include <random>

using namespace ausaxs;
using namespace ausaxs::rigidbody::parameter;

std::tuple<
    std::uniform_real_distribution<double>, std::uniform_real_distribution<double>, 
    std::uniform_real_distribution<double>, std::uniform_real_distribution<double>
> create_distributions(double length_start, double rad_start, double Rg) {
    std::uniform_real_distribution<double> translation_dist(-length_start, length_start);
    std::uniform_real_distribution<double> rotation_dist(-rad_start, rad_start);
    std::uniform_real_distribution<double> translation_symmetry_dist(-Rg, Rg);
    std::uniform_real_distribution<double> rotation_symmetry_dist(-3, 3);
    return {std::move(translation_dist), std::move(rotation_dist), std::move(translation_symmetry_dist), std::move(rotation_symmetry_dist)};
}

ParameterGenerationStrategy::ParameterGenerationStrategy(
    observer_ptr<const Rigidbody> rigidbody, unsigned int iterations, double length_start, double rad_start) 
    : rigidbody(rigidbody), decay_strategy(rigidbody::factory::create_decay_strategy(iterations)
) {
    std::tie(translation_dist, rotation_dist, translation_symmetry_dist, rotation_symmetry_dist) = create_distributions(length_start, rad_start, rigidbody->molecule.get_Rg());
}

ParameterGenerationStrategy::ParameterGenerationStrategy(
    observer_ptr<const Rigidbody> rigidbody, std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy, double length_start, double rad_start) 
    : rigidbody(rigidbody), decay_strategy(std::move(decay_strategy)
) {
    std::tie(translation_dist, rotation_dist, translation_symmetry_dist, rotation_symmetry_dist) = create_distributions(length_start, rad_start, rigidbody->molecule.get_Rg());
}

ParameterGenerationStrategy::~ParameterGenerationStrategy() = default;

void ParameterGenerationStrategy::set_decay_strategy(std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy) {
    this->decay_strategy = std::move(decay_strategy);
}

observer_ptr<rigidbody::parameter::decay::DecayStrategy> ParameterGenerationStrategy::get_decay_strategy() const {return decay_strategy.get();}

void ParameterGenerationStrategy::set_max_translation_distance(double distance) {
    translation_dist = std::uniform_real_distribution<double>(-distance, distance);
}

void ParameterGenerationStrategy::set_max_rotation_angle(double radians) {
    rotation_dist = std::uniform_real_distribution<double>(-radians, radians);
}