// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/CopyBodyElement.h>
#include <rigidbody/sequencer/detail/BodyIndexOps.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/selection/SymmetryTargets.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/observer_ptr.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

void clone(observer_ptr<Sequencer> owner, std::string_view body_name, int index) {
    ausaxs::rigidbody::sequencer::detail::require_mutable_structure(owner, "copy");

    // initial_conformation stores bodies centered at origin; the live molecule body stores absolute positions
    auto initial_body = owner->_get_rigidbody()->conformation->initial_conformation[index];
    auto body_pars = owner->_get_rigidbody()->conformation->absolute_parameters.parameters[index];
    body_pars.translation.z() += 2*owner->_get_molecule()->get_Rg(false);

    owner->_get_molecule()->get_bodies().emplace_back(owner->_get_molecule()->get_body(index));
    owner->_get_rigidbody()->conformation->initial_conformation.emplace_back(std::move(initial_body));
    owner->_get_rigidbody()->conformation->absolute_parameters.parameters.emplace_back(body_pars);

    int new_index = static_cast<int>(owner->_get_molecule()->size_body())-1;
    owner->setup()._body_name_registry().add_body(new_index, std::string{body_name});
    owner->_get_rigidbody()->symmetry_targets->invalidate(); // the copy brings its source's symmetries with it
    owner->_get_rigidbody()->constraints->invalidate();      // the new body needs an entry of its own in the per-body constraint map
}

CopyBodyElement::CopyBodyElement(observer_ptr<Sequencer> owner, std::string_view body_name, std::string_view source_body_name) {
    clone(owner, body_name, owner->setup()._get_body(source_body_name));
}

CopyBodyElement::CopyBodyElement(observer_ptr<Sequencer> owner, std::string_view body_name, int source_body_index) {
    clone(owner, body_name, source_body_index);
}

CopyBodyElement::~CopyBodyElement() = default;

void CopyBodyElement::run() {}

std::vector<std::string> CopyBodyElement::_valid_arguments() {
    return {};
}

std::unique_ptr<GenericElement> CopyBodyElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (args.inlined.size() != 2) {throw except::parse_error(
        "copy", "Invalid number of inline arguments. Expected [new name] [target name], but got " + std::to_string(args.inlined.size()) + "."
    );}

    const auto& body_names = owner->_get_sequencer()->setup()._body_name_registry();
    std::string source = args.inlined[0];
    std::string name = args.inlined[1];
    if (!body_names.contains(source)) {
        if (body_names.contains(name)) {std::swap(source, name);}
        else {throw except::parse_error("copy", "Body name \"" + source + "\" not found.");}
    }
    if (body_names.contains(name)) {throw except::parse_error("copy", "Body name \"" + name + "\" already exists.");}

    return std::make_unique<CopyBodyElement>(
        owner->_get_sequencer(),
        name,
        source
    );
}