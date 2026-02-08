// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/SetupElement.h>
#include <rigidbody/sequencer/elements/setup/LoadExistingElement.h>
#include <rigidbody/sequencer/elements/setup/LoadElement.h>
#include <rigidbody/sequencer/elements/setup/AutoConstraintsElement.h>
#include <rigidbody/sequencer/elements/setup/SymmetryElement.h>
#include <rigidbody/sequencer/elements/setup/RelativeHydrationElement.h>
#include <rigidbody/sequencer/elements/setup/OutputFolderElement.h>
#include <rigidbody/sequencer/elements/setup/ConstraintElement.h>
#include <rigidbody/sequencer/elements/LoopElementCallback.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/Rigidbody.h>

using namespace ausaxs::rigidbody::sequencer;

SetupElement::SetupElement(observer_ptr<Sequencer> owner) : LoopElementCallback(owner) {}
SetupElement::SetupElement(observer_ptr<Sequencer> owner, io::ExistingFile saxs) : LoopElementCallback(owner), saxs_path(std::move(saxs)) {}

LoopElement& SetupElement::end() {
    return *owner;
}

SetupElement& SetupElement::set_overlap_function(std::function<double(double)> func) {
    rigidbody::constraints::OverlapConstraint::set_overlap_function(std::move(func));
    return *this;
}

SetupElement& SetupElement::load(const std::vector<std::string>& paths, const std::vector<std::string>& names) {
    elements.push_back(std::make_unique<LoadElement>(owner->_get_sequencer(), paths, names));
    return *this;
}

SetupElement& SetupElement::load(const std::string& path, const std::vector<int>& split, const std::vector<std::string>& names) {
    elements.push_back(std::make_unique<LoadElement>(owner->_get_sequencer(), path, split, names));
    return *this;
}

SetupElement& SetupElement::load(const std::string& path, const std::vector<std::string>& names) {
    elements.push_back(std::make_unique<LoadElement>(owner->_get_sequencer(), path, names));
    return *this;
}

SetupElement& SetupElement::load_existing(observer_ptr<Rigidbody> rigidbody) {
    elements.push_back(std::make_unique<LoadExistingElement>(owner->_get_sequencer(), rigidbody));
    return *this;
}

std::unordered_map<std::string, unsigned int>& SetupElement::_get_body_names() {
    return body_names;
}

detail::BodySymmetrySelector SetupElement::_get_body_index(std::string_view name) const {
    auto it = body_names.find(std::string{name});
    if (it == body_names.end()) {
        throw std::runtime_error(
            "SetupElement::_get_body_index: Unknown body name \"" + std::string{name} + "\". "
            "Known body names are: " + [&] () {
                std::string result;
                for (const auto& [name, index] : body_names) {
                    result += name + " ";
                }
                return result;
            }()
        );
    }
    return detail::from_index(it->second);
}

void SetupElement::_set_active_body(observer_ptr<Rigidbody> body) {
    active_body = body;
    owner->_get_sequencer()->rigidbody = body;
}

SetupElement& SetupElement::distance_constraint(const std::string& body1, const std::string& body2, unsigned int iatom1, unsigned int iatom2) {
    owner->_get_rigidbody()->constraints->add_constraint(
        std::make_unique<constraints::DistanceConstraintAtom>(
            &active_body->molecule,
            body_names.at(body1),
            body_names.at(body2),
            iatom1,
            iatom2
        )
    );
    return *this;
}

SetupElement& SetupElement::distance_constraint_closest(unsigned int ibody1, unsigned int ibody2) {
    owner->_get_rigidbody()->constraints->add_constraint(
        std::make_unique<constraints::DistanceConstraintBond>(
            &active_body->molecule,
            ibody1, 
            ibody2
        )
    );
    return *this;
}

SetupElement& SetupElement::distance_constraint_closest(const std::string& ibody1, const std::string& ibody2) {
    owner->_get_rigidbody()->constraints->add_constraint(
        std::make_unique<constraints::DistanceConstraintBond>(
            &active_body->molecule,
            body_names.at(ibody1), 
            body_names.at(ibody2)
        )
    );
    return *this;
}

SetupElement& SetupElement::distance_constraint_center_mass(unsigned int ibody1, unsigned int ibody2) {
    owner->_get_rigidbody()->constraints->add_constraint(
        std::make_unique<constraints::DistanceConstraintCM>(
            &active_body->molecule,
            ibody1, 
            ibody2
        )
    );
    return *this;
}

SetupElement& SetupElement::distance_constraint_center_mass(const std::string& ibody1, const std::string& ibody2) {
    owner->_get_rigidbody()->constraints->add_constraint(
        std::make_unique<constraints::DistanceConstraintCM>(
            &active_body->molecule,
            body_names.at(ibody1), 
            body_names.at(ibody2)
        )
    );
    return *this;
}

std::string SetupElement::_get_config_folder() const {
    return config_folder;
}

void SetupElement::_set_config_folder(const io::Folder& folder) {
    config_folder = folder;
}

void SetupElement::_set_saxs_path(const io::ExistingFile& saxs) {
    saxs_path = saxs;
}

SetupElement& SetupElement::fixed_constraint() {
    throw std::runtime_error("SetupElement::fixed_constraint: Not implemented.");
}

SetupElement& SetupElement::generate_linear_constraints() {
    elements.push_back(std::make_unique<AutoConstraintsElement>(static_cast<Sequencer*>(owner), settings::rigidbody::ConstraintGenerationStrategyChoice::Linear));
    return *this;
}

SetupElement& SetupElement::generate_volumetric_constraints() {
    elements.push_back(std::make_unique<AutoConstraintsElement>(static_cast<Sequencer*>(owner), settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric));
    return *this;
}

SetupElement& SetupElement::add_constraint(std::unique_ptr<rigidbody::constraints::Constraint> constraint) {
    elements.push_back(std::make_unique<ConstraintElement>(owner->_get_sequencer(), std::move(constraint)));
    return *this;
}

SetupElement& SetupElement::symmetry(const std::vector<std::string>& names, const std::vector<symmetry::type>& symmetry) {
    elements.push_back(std::make_unique<SymmetryElement>(owner->_get_sequencer(), names, symmetry));
    return *this;
}

SetupElement& SetupElement::symmetry(std::string_view name, symmetry::type symmetry) {
    this->symmetry(std::vector<std::string>{std::string{name}}, std::vector<symmetry::type>{symmetry});
    return *this;
}

SetupElement& SetupElement::relative_hydration(const std::vector<std::string>& names, const std::vector<double>& ratios) {
    elements.push_back(std::make_unique<RelativeHydrationElement>(owner->_get_sequencer(), names, ratios));
    return *this;
}

SetupElement& SetupElement::output_folder(const io::Folder& path) {
    elements.push_back(std::make_unique<OutputFolderElement>(owner->_get_sequencer(), path));
    return *this;
}

const ausaxs::io::ExistingFile& SetupElement::_get_saxs_path() const {
    return saxs_path;
}

std::vector<std::unique_ptr<GenericElement>>& SetupElement::_get_elements() {
    return elements;
}