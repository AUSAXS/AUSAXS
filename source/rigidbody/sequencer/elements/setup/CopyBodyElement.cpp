// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/CopyBodyElement.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/observer_ptr.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

void clone(observer_ptr<Sequencer> owner, std::string_view body_name, int index) {
    // initial_conformation stores bodies centered at origin; the live molecule body stores absolute positions
    auto initial_body = owner->_get_rigidbody()->conformation->initial_conformation[index];
    auto body_pars = owner->_get_rigidbody()->conformation->absolute_parameters.parameters[index];
    body_pars.translation.z() += 2*owner->_get_molecule()->get_Rg(false);

    owner->_get_molecule()->get_bodies().emplace_back(owner->_get_molecule()->get_body(index));
    owner->_get_rigidbody()->conformation->initial_conformation.emplace_back(std::move(initial_body));
    owner->_get_rigidbody()->conformation->absolute_parameters.parameters.emplace_back(body_pars);

    unsigned int new_index = owner->_get_molecule()->size_body()-1;
    owner->setup()._get_body_names().emplace(std::string{body_name}, rigidbody::sequencer::detail::to_index(new_index));
}

CopyBodyElement::CopyBodyElement(observer_ptr<Sequencer> owner, std::string_view body_name, std::string_view source_body_name) {
    auto index = owner->setup()._get_body_index(source_body_name);
    clone(owner, body_name, index.body);
}

CopyBodyElement::CopyBodyElement(observer_ptr<Sequencer> owner, std::string_view body_name, int source_body_index) {
    clone(owner, body_name, source_body_index);
}

CopyBodyElement::~CopyBodyElement() = default;

void CopyBodyElement::run() {}