// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/SetupElement.h>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/sequencer/Sequencer.h>

using namespace ausaxs::rigidbody::sequencer;

SetupElement::SetupElement(observer_ptr<Sequencer> sequencer) : owner(sequencer) {}
SetupElement::SetupElement(observer_ptr<Sequencer> sequencer, io::ExistingFile saxs) : owner(sequencer), saxs_path(std::move(saxs)) {}

detail::BodyNameRegistry& SetupElement::_body_name_registry() {
    return body_names;
}

detail::BodySymmetrySelector SetupElement::_get_body_index(std::string_view name) const {
    return body_names.resolve(name);
}

int SetupElement::_get_body(std::string_view name) const {
    return body_names.resolve_body(name);
}

void SetupElement::_set_active_body(observer_ptr<Rigidbody> body) {
    active_body = body;
    owner->_get_sequencer()->rigidbody = body;
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

const ausaxs::io::ExistingFile& SetupElement::_get_saxs_path() const {
    return saxs_path;
}

std::vector<std::unique_ptr<GenericElement>>& SetupElement::_get_elements() {
    return elements;
}