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
#include <data/Water.h>

#include <iostream>

using namespace rigidbody::sequencer;

template<typename T> requires std::is_same_v<std::decay_t<T>, Protein>
RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, T&& rigidbody) : RigidBody(std::forward<T>(rigidbody)) {
    prepare_fitter(saxs);
    initialize();
}
template RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, Protein&& rigidbody);
template RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, const Protein& rigidbody);
RigidBodyManager::~RigidBodyManager() = default;

void RigidBodyManager::initialize() {
    rigidbody->generate_new_hydration();

    // save the best configuration in a simple struct
    best = detail::BestConf {
        .grid = std::make_shared<grid::Grid>(*rigidbody->get_grid()),
        .waters = waters(),
        .chi2 = fitter->fit_only()
    };
}

void RigidBodyManager::optimize_step() {
    RigidBody::optimize_step(best);
}

void RigidBodyManager::set_managers(const settings::rigidbody::BodySelectStrategyChoice& body_selector, const settings::rigidbody::TransformationStrategyChoice& transform, const settings::rigidbody::ParameterGenerationStrategyChoice& parameters) {
    this->body_selector = std::move(rigidbody::factory::create_selection_strategy(rigidbody.get(), body_selector));
    this->transform = std::move(rigidbody::factory::create_transform_strategy(rigidbody.get(), transform));
    this->parameter_generator = std::move(rigidbody::factory::create_parameter_strategy(settings::rigidbody::iterations, 5, M_PI/3, parameters));
}