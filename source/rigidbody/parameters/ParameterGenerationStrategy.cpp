/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <math/Vector3.h>

#include <random>

using namespace ausaxs::rigidbody::parameter;

ParameterGenerationStrategy::ParameterGenerationStrategy(
    observer_ptr<const RigidBody> molecule, unsigned int iterations, double length_start, double rad_start) 
    : molecule(molecule), decay_strategy(rigidbody::factory::create_decay_strategy(iterations)
) {
    generator = std::mt19937(std::random_device{}());
    translation_dist = std::uniform_real_distribution<double>(-length_start, length_start);
    rotation_dist = std::uniform_real_distribution<double>(-rad_start, rad_start);
    symmetry_dist = std::uniform_real_distribution<double>(-10, 10);
}

ParameterGenerationStrategy::ParameterGenerationStrategy(
    observer_ptr<const RigidBody> molecule, std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy, double length_start, double rad_start) 
    : molecule(molecule), decay_strategy(std::move(decay_strategy)
) {
    generator = std::mt19937(std::random_device{}());
    translation_dist = std::uniform_real_distribution<double>(-length_start, length_start);
    rotation_dist = std::uniform_real_distribution<double>(-rad_start, rad_start);
    symmetry_dist = std::uniform_real_distribution<double>(-10, 10);
}

ParameterGenerationStrategy::~ParameterGenerationStrategy() = default;

void ParameterGenerationStrategy::set_decay_strategy(std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy) {
    this->decay_strategy = std::move(decay_strategy);
}

void ParameterGenerationStrategy::set_max_translation_distance(double distance) {
    translation_dist = std::uniform_real_distribution<double>(-distance, distance);
}

void ParameterGenerationStrategy::set_max_rotation_angle(double radians) {
    rotation_dist = std::uniform_real_distribution<double>(-radians, radians);
}