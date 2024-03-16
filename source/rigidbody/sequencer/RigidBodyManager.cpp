/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/RigidBodyManager.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <hydrate/Grid.h>
#include <hydrate/GridMember.h>
#include <hydrate/placement/PlacementStrategy.h>
#include <hydrate/culling/CullingStrategy.h>
#include <io/ExistingFile.h>
#include <fitter/LinearFitter.h>
#include <data/record/Water.h>
#include <constants/ConstantsMath.h>

#include <iostream>

using namespace rigidbody::sequencer;

std::unique_ptr<RigidBodyManager> rigidbody::sequencer::rigidbody;

template<typename T> requires std::is_same_v<std::decay_t<T>, data::Molecule>
RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, T&& rigidbody) : RigidBody(std::forward<T>(rigidbody)) {
    prepare_fitter(saxs);
    initialize();
}
template RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, data::Molecule&& rigidbody);
template RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, const data::Molecule& rigidbody);
RigidBodyManager::~RigidBodyManager() = default;

void RigidBodyManager::initialize() {
    rigidbody->generate_new_hydration();

    // save the best configuration in a simple struct
    best = detail::BestConf(std::make_shared<grid::Grid>(*rigidbody->get_grid()), get_waters(), fitter->fit_chi2_only());
}

void RigidBodyManager::optimize_step() {
    RigidBody::optimize_step(best);
}

void RigidBodyManager::set_managers(const settings::rigidbody::BodySelectStrategyChoice& body_selector, const settings::rigidbody::TransformationStrategyChoice& transform, const settings::rigidbody::ParameterGenerationStrategyChoice& parameters) {
    this->body_selector = rigidbody::factory::create_selection_strategy(rigidbody.get(), body_selector);
    this->transform = rigidbody::factory::create_transform_strategy(rigidbody.get(), transform);
    this->parameter_generator = rigidbody::factory::create_parameter_strategy(settings::rigidbody::iterations, 5, constants::pi/3, parameters);
}