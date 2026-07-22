// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/DeleteElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/detail/BodyIndexOps.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <utility/observer_ptr.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

DeleteElement::DeleteElement(observer_ptr<Sequencer> owner, std::vector<std::string> names) {
    std::vector<int> indices;
    indices.reserve(names.size());
    for (const auto& name : names) {
        indices.push_back(owner->setup()._get_body(name));
    }
    detail::erase_bodies(owner, std::move(indices));
}

DeleteElement::~DeleteElement() = default;

void DeleteElement::run() {}

std::vector<std::string> DeleteElement::_valid_arguments() {
    return {};
}

std::unique_ptr<GenericElement> DeleteElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("delete", "Unexpected named argument \"" + args.named.begin()->first + "\".");}
    if (args.inlined.empty()) {throw except::parse_error(
        "delete", "Invalid number of inline arguments. Expected one or more body names, but got 0."
    );}

    const auto& body_names = owner->_get_sequencer()->setup()._get_body_names();
    std::vector<std::string> names;
    names.reserve(args.inlined.size());
    for (std::size_t i = 0; i < args.inlined.size(); ++i) {
        const std::string& name = args.inlined[i];
        if (!body_names.contains(name)) {throw except::parse_error("delete", "Body name \"" + name + "\" not found.");}
        if (std::find(names.begin(), names.end(), name) != names.end()) {
            throw except::parse_error("delete", "Body name \"" + name + "\" was specified more than once.");
        }
        names.push_back(name);
    }

    auto total_bodies = owner->_get_sequencer()->_get_molecule()->size_body();
    if (names.size() >= total_bodies) {throw except::parse_error("delete", "Cannot delete all bodies; at least one must remain.");}

    return std::make_unique<DeleteElement>(owner->_get_sequencer(), std::move(names));
}
