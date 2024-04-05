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

using namespace rigidbody::sequencer;

std::unique_ptr<RigidBodyManager> rigidbody::sequencer::rigidbody;

template<typename T> requires std::is_base_of_v<data::Molecule, std::decay_t<T>>
RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, T&& rigidbody) : RigidBody(std::forward<T>(rigidbody)) {
    initialize();
    prepare_fitter(saxs);
}
RigidBodyManager::~RigidBodyManager() = default;

void RigidBodyManager::set_constraint_manager(std::shared_ptr<rigidbody::constraints::ConstraintManager> constraints) {
    this->constraints = constraints;
}

void RigidBodyManager::set_body_select_manager(std::shared_ptr<rigidbody::selection::BodySelectStrategy> body_selector) {
    this->body_selector = body_selector;
}

void RigidBodyManager::set_transform_manager(std::shared_ptr<rigidbody::transform::TransformStrategy> transform) {
    this->transform = transform;
}

void RigidBodyManager::set_parameter_manager(std::shared_ptr<rigidbody::parameter::ParameterGenerationStrategy> parameters) {
    this->parameter_generator = parameters;
}

void RigidBodyManager::initialize() {
    rigidbody->generate_new_hydration();
    best = detail::BestConf(std::make_shared<grid::Grid>(*rigidbody->get_grid()), get_waters(), fitter->fit_chi2_only());
}

void RigidBodyManager::optimize_step() {
    RigidBody::optimize_step(best);
}

std::shared_ptr<fitter::Fit> RigidBodyManager::get_fit() const {
    return fitter->get_fit();
}

template RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, data::Molecule&& rigidbody);
template RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, const data::Molecule& rigidbody);
template RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, data::Molecule& rigidbody);
template RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, rigidbody::RigidBody&& rigidbody);
template RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, const rigidbody::RigidBody& rigidbody);
template RigidBodyManager::RigidBodyManager(const io::ExistingFile& saxs, rigidbody::RigidBody& rigidbody);