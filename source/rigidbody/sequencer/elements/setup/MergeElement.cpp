// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/MergeElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/observer_ptr.h>

#include <algorithm>
#include <iterator>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    // dst is kept only if both sides track it; otherwise the merged result can no longer be trusted to be parallel-indexed with the atoms
    template<typename T>
    void merge_optional_vector(std::optional<std::vector<T>>& dst, const std::optional<std::vector<T>>& src) {
        if (dst.has_value() && src.has_value()) {dst->insert(dst->end(), src->begin(), src->end());}
        else {throw std::runtime_error("MergeElement: metadata mismatch. This is not supposed to happen.");}
    }

    void merge_metadata(data::Body& first, const data::Body& other) {
        if (!first.get_metadata().has_value() && !other.get_metadata().has_value()) {return;}
        data::AtomMetadata result = first.get_metadata().value_or(data::AtomMetadata{});
        data::AtomMetadata other_meta = other.get_metadata().value_or(data::AtomMetadata{});
        merge_optional_vector(result.backbone, other_meta.backbone);
        merge_optional_vector(result.residue_seq, other_meta.residue_seq);
        merge_optional_vector(result.occupancy, other_meta.occupancy);
        first.set_metadata(std::move(result));
    }
}

MergeElement::MergeElement(observer_ptr<Sequencer> owner, std::string_view first_name, std::vector<std::string> other_names) {
    auto& body_names = owner->setup()._get_body_names();
    int i_first = owner->setup()._get_body_index(first_name).body;

    std::vector<int> other_indices;
    other_indices.reserve(other_names.size());
    for (const auto& name : other_names) {
        other_indices.push_back(owner->setup()._get_body_index(name).body);
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

    // remove the merged-away bodies, highest index first so earlier indices stay valid while erasing
    std::vector<int> removed = other_indices;
    std::sort(removed.begin(), removed.end());
    for (auto it = removed.rbegin(); it != removed.rend(); ++it) {
        molecule.get_bodies().erase(molecule.get_bodies().begin() + *it);
        conformation.initial_conformation.erase(conformation.initial_conformation.begin() + *it);
        conformation.absolute_parameters.parameters.erase(conformation.absolute_parameters.parameters.begin() + *it);
    }

    // drop the merged-away names and reindex the survivors to their shifted body indices
    for (auto it = body_names.begin(); it != body_names.end(); ) {
        auto sel = detail::from_index(it->second);
        auto pos = std::lower_bound(removed.begin(), removed.end(), sel.body);
        if (pos != removed.end() && *pos == sel.body) {
            it = body_names.erase(it);
            continue;
        }
        auto shift = std::distance(removed.begin(), pos);
        if (shift != 0) {it->second = detail::to_index(sel.body - static_cast<int>(shift), sel.symmetry, sel.replica);}
        ++it;
    }
}

MergeElement::~MergeElement() = default;

void MergeElement::run() {}

std::vector<std::string> MergeElement::_valid_arguments() {
    return {};
}

std::unique_ptr<GenericElement> MergeElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("merge", "Unexpected named argument \"" + args.named.begin()->first + "\".");}
    if (args.inlined.size() < 2) {throw except::parse_error(
        "merge", "Invalid number of inline arguments. Expected [first] [others...], but got " + std::to_string(args.inlined.size()) + "."
    );}

    const auto& body_names = owner->_get_sequencer()->setup()._get_body_names();
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
