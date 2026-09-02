// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/MergeElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/detail/BodyIndexOps.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/observer_ptr.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    void merge_metadata(data::Body& first, const data::Body& other) {
        if (!first.get_metadata().has_value() && !other.get_metadata().has_value()) {return;}
        data::AtomMetadata result = first.get_metadata().value_or(data::AtomMetadata{});

        // append rejects a field engaged on only one side, which would mean the two bodies were tracked under different settings. That should not be reachable
        // from the sequencer, as the setting is fixed for the duration of a run.
        result.append(other.get_metadata().value_or(data::AtomMetadata{}));
        first.set_metadata(std::move(result));
    }
}

MergeElement::MergeElement(observer_ptr<Sequencer> owner, std::string_view first_name, std::vector<std::string> other_names) {
    detail::require_mutable_structure(owner, "merge");
    int i_first = owner->setup()._get_body(first_name);

    std::vector<int> other_indices;
    other_indices.reserve(other_names.size());
    for (const auto& name : other_names) {
        other_indices.push_back(owner->setup()._get_body(name));
    }

    auto& molecule = *owner->_get_molecule();
    auto& conformation = *owner->_get_rigidbody()->conformation;
    auto& first_body = molecule.get_body(i_first);

    for (int idx : other_indices) {
        auto& other_body = molecule.get_body(idx);
        merge_metadata(first_body, other_body);

        auto& first_atoms = first_body.get_atoms();
        const auto& other_atoms = other_body.get_atoms();
        first_atoms.insert(first_atoms.end(), other_atoms.begin(), other_atoms.end());
    }

    // recompute the merged body's center of mass, refreshing the origin-centered initial conformation to match
    auto cm = first_body.get_cm();
    data::Body centered = first_body;
    centered.translate(-cm);
    conformation.initial_conformation[i_first] = std::move(centered);
    conformation.absolute_parameters.parameters[i_first].translation = cm;

    detail::erase_bodies(owner, std::move(other_indices));

    // the body count changed, so the histogram manager's per-body-indexed caches (sized for the old body count) must be rebuilt from scratch. The grid
    // is discarded rather than refreshed: Molecule::get_grid rebuilds it on demand, so doing it here would only force the work at parse time.
    owner->_get_molecule()->reset_histogram_manager();
    owner->_get_rigidbody()->molecule.clear_grid();
}

MergeElement::~MergeElement() = default;

void MergeElement::run() {}

std::vector<std::string> MergeElement::_valid_arguments() {
    return {};
}

InlineSignature MergeElement::_valid_inline_arguments() {
    return {.names = {"first", "others..."}, .min = 2, .max = unbounded_inline_args};
}

// merge [first] [others...] - merges every [others] body into [first]
std::unique_ptr<GenericElement> MergeElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    const auto& body_names = owner->_get_sequencer()->setup()._body_name_registry();
    std::string first = args.inlined[0];
    if (!body_names.contains(first)) {throw except::parse_error("merge", "Body name \"" + first + "\" not found.");}

    std::vector<std::string> others;
    others.reserve(args.inlined.size() - 1);
    for (std::size_t i = 1; i < args.inlined.size(); ++i) {
        const std::string& name = args.inlined[i];
        if (name == first) {throw except::parse_error("merge", "Cannot merge body \"" + name + "\" into itself.");}
        if (!body_names.contains(name)) {throw except::parse_error("merge", "Body name \"" + name + "\" not found.");}
        if (std::find(others.begin(), others.end(), name) != others.end()) {
            throw except::parse_error("merge", "Body name \"" + name + "\" was specified more than once.");
        }
        others.push_back(name);
    }

    return std::make_unique<MergeElement>(
        owner->_get_sequencer(),
        first,
        std::move(others)
    );
}
