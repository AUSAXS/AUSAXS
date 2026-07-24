// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/SplitElement.h>
#include <rigidbody/sequencer/detail/BodyIndexOps.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <rigidbody/Rigidbody.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <data/atoms/AtomMetadata.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/StringUtils.h>

#include <algorithm>
#include <cassert>
#include <limits>
#include <optional>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    // Partition `body`'s atoms (and, in lockstep, its metadata) into splits.size()+1 fragments at the given residue
    // sequence ids. Mirrors BodySplitter::split's boundary logic, but reads residue ids from the body's own
    // per-atom metadata rather than re-reading a source file - the body may already have been transformed, merged,
    // or had its symmetry converted since it was loaded, so there is no file to consistently re-split.
    //
    // Fragments are built from atoms only: any existing explicit hydration (e.g. crystallographic waters loaded
    // from the source file) on `body` has no well-defined per-atom-range assignment across fragments, so it is
    // dropped - matching BodySplitter's own load-time split, which drops waters at split boundaries the same way.
    std::vector<data::Body> partition_by_residue(const data::Body& body, const std::vector<int>& splits) {
        const auto& metadata = body.get_metadata();
        if (!metadata || !metadata->residue_seq) {
            throw ausaxs::rigidbody::sequencer::except::parse_error("split",
                "Body has no residue sequence metadata; enable settings::molecule::store_residue_seq before "
                "loading to make splitting by residue possible (the sequencer already enables this by default)."
            );
        }
        const auto& atoms = body.get_atoms();
        const auto& resseq = *metadata->residue_seq;
        assert(resseq.size() == atoms.size() && "partition_by_residue: residue_seq metadata is not parallel-indexed to the atom vector.");

        int max_id = std::numeric_limits<int>::min(), min_id = std::numeric_limits<int>::max();
        for (int r : resseq) {
            min_id = std::min(min_id, r);
            max_id = std::max(max_id, r);
        }

        // 1 extra to allow splitting after the last residue, another for 0-based indexing
        std::vector<bool> split_at(max_id+2, false);
        for (int id : splits) {
            if (id < 0 || id < min_id) {
                throw ausaxs::rigidbody::sequencer::except::parse_error("split", "Split " + std::to_string(id) + " smaller than lowest residue sequence id (" + std::to_string(min_id) + ").");
            }
            if (id > max_id+1) {
                throw ausaxs::rigidbody::sequencer::except::parse_error("split", "Split " + std::to_string(id) + " larger than highest residue sequence id (" + std::to_string(max_id) + ").");
            }
            split_at[id] = true;
        }

        auto slice_metadata = [&metadata] (std::size_t begin, std::size_t end) -> std::optional<data::AtomMetadata> {
            data::AtomMetadata m;
            bool any = false;
            if (metadata->backbone)    {m.backbone    = std::vector<data::backbone_t>(metadata->backbone->begin()+begin, metadata->backbone->begin()+end); any = true;}
            if (metadata->residue_seq) {m.residue_seq = std::vector<int>(metadata->residue_seq->begin()+begin, metadata->residue_seq->begin()+end);        any = true;}
            if (metadata->occupancy)   {m.occupancy   = std::vector<float>(metadata->occupancy->begin()+begin, metadata->occupancy->begin()+end);          any = true;}
            return any ? std::make_optional(std::move(m)) : std::nullopt;
        };

        std::vector<data::Body> fragments(splits.size()+1);
        std::size_t index_body = 0;
        std::size_t begin = 0;
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            int r = std::max(resseq[i], 0); // in some files resSeq starts negative
            if (split_at[r]) {
                fragments[index_body] = data::Body(std::vector<data::AtomFF>(atoms.begin()+begin, atoms.begin()+i));
                if (auto m = slice_metadata(begin, i)) {fragments[index_body].set_metadata(std::move(*m));}
                ++index_body;
                split_at[r] = false; // mark it as false so we won't split again on the next atom
                begin = i;
            }
        }
        fragments[index_body] = data::Body(std::vector<data::AtomFF>(atoms.begin()+begin, atoms.end()));
        if (auto m = slice_metadata(begin, atoms.size())) {fragments[index_body].set_metadata(std::move(*m));}
        return fragments;
    }

    // A freshly-constructed Body carries a plain (non-optimizable) SymmetryStorage; every body already in a
    // Rigidbody was converted once, at construction time (see Rigidbody::Rigidbody). New fragments need the
    // same conversion before any symmetry can be attached and driven by the optimiser.
    void make_optimizable(data::Body& body) {
        if (dynamic_cast<symmetry::OptimizableSymmetryStorage*>(body.symmetry().get_obj()) == nullptr) {
            body.symmetry().set_obj(std::make_unique<symmetry::OptimizableSymmetryStorage>(std::move(*body.symmetry().get_obj())));
        }
    }
}

SplitElement::SplitElement(observer_ptr<Sequencer> owner, const std::string& body_name, std::vector<int> splits) : owner(owner) {
    _split(body_name, splits);
}

SplitElement::~SplitElement() = default;

void SplitElement::run() {}

void SplitElement::_split(const std::string& body_name, const std::vector<int>& splits) {
    auto molecule = owner->_get_molecule();
    auto rigidbody = owner->_get_rigidbody();
    auto& setup = owner->setup();

    auto index = setup._get_body_index(body_name);
    if (index.symmetry != -1 || index.replica != 0) {
        throw except::parse_error("split", "Body \"" + body_name + "\" is a symmetry replica, not a base body.");
    }
    int ib = index.body;
    const auto& original = molecule->get_body(ib);

    // splitting a body that already participates in a symmetry shared with bodies outside the split would require
    // expanding that ReferenceSymmetry's membership (or reassigning its ownership, if ib is the primary) - not yet
    // supported, so reject rather than silently produce an inconsistent shared symmetry
    for (std::size_t i = 0; i < original.size_symmetry(); ++i) {
        if (dynamic_cast<const symmetry::ReferenceSymmetryView*>(original.symmetry().get(i))) {
            throw except::parse_error("split", "Body \"" + body_name + "\" participates in a symmetry shared with other bodies; splitting it is not yet supported.");
        }
    }
    for (std::size_t b = 0; b < molecule->size_body(); ++b) {
        if (static_cast<int>(b) == ib) {continue;}
        const auto& other = molecule->get_body(b);
        for (std::size_t i = 0; i < other.size_symmetry(); ++i) {
            auto* view = dynamic_cast<const symmetry::ReferenceSymmetryView*>(other.symmetry().get(i));
            if (view != nullptr && view->primary_body == ib) {
                throw except::parse_error("split", "Body \"" + body_name + "\" owns a symmetry shared with other bodies; splitting it is not yet supported.");
            }
        }
    }

    // clone the body's existing symmetries (if any) before partitioning it away; each becomes a ReferenceSymmetry
    // shared by every resulting fragment, rather than an independently-optimizable copy on each
    std::vector<std::unique_ptr<symmetry::ISymmetry>> orig_symmetries;
    orig_symmetries.reserve(original.size_symmetry());
    for (std::size_t i = 0; i < original.size_symmetry(); ++i) {
        orig_symmetries.push_back(original.symmetry().get(i)->clone());
    }

    auto fragments = partition_by_residue(original, splits);
    if (fragments.size() < 2) {throw except::parse_error("split", "Expected at least one split index.");}
    for (auto& frag : fragments) {make_optimizable(frag);}

    // append the fragments at the tail of the molecule/conformation; erasing the original body afterwards shifts
    // every one of these indices down by exactly one, since they are necessarily all above it
    std::size_t append_at = molecule->size_body();
    for (auto& frag : fragments) {
        molecule->get_bodies().emplace_back(std::move(frag));
        rigidbody->conformation->initial_conformation.emplace_back(); // placeholder; overwritten once symmetries are installed
        rigidbody->conformation->absolute_parameters.parameters.emplace_back();
    }
    detail::erase_bodies(owner, {ib});

    std::vector<int> final_indices(fragments.size());
    for (std::size_t k = 0; k < fragments.size(); ++k) {final_indices[k] = static_cast<int>(append_at + k) - 1;}

    auto& name_map = setup._body_name_registry();
    for (int idx : final_indices) {name_map.add_body(idx);}

    // install each of the original symmetries as a ReferenceSymmetry shared by every fragment: the first fragment
    // becomes the owning primary, the rest hold non-owning views (same pattern as SymmetryElement::_add_reference)
    int primary = final_indices.front();
    for (std::size_t i = 0; i < orig_symmetries.size(); ++i) {
        auto* opt_primary = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(molecule->get_body(primary).symmetry().get_obj());
        opt_primary->optimize_translate = true;
        opt_primary->optimize_rot_axis = true;
        molecule->get_body(primary).symmetry().add(std::make_unique<symmetry::ReferenceSymmetry>(
            orig_symmetries[i]->clone(), final_indices, std::vector<int>(final_indices.size(), static_cast<int>(i)), molecule
        ));

        auto* ref = static_cast<symmetry::ReferenceSymmetry*>(molecule->get_body(primary).symmetry().get(i));
        int reps = static_cast<int>(ref->repetitions());

        for (std::size_t k = 1; k < final_indices.size(); ++k) {
            int b = final_indices[k];
            auto* opt = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(molecule->get_body(b).symmetry().get_obj());
            opt->optimize_translate = true;
            opt->optimize_rot_axis = true;
            molecule->get_body(b).symmetry().add(std::make_unique<symmetry::ReferenceSymmetryView>(molecule, primary, static_cast<int>(i)));
        }

        for (int b : final_indices) {
            for (int j = 0; j < reps; ++j) {name_map.add_replica(b, static_cast<int>(i), j+1);}
        }
    }

    // derive each fragment's initial (origin-centred) conformation fresh from its current, now symmetry-equipped
    // live state - mirrors SystemSpecification's own bootstrap recipe, since these are brand-new bodies with no
    // prior initial_conformation entry that needs to stay in sync
    for (int idx : final_indices) {
        data::Body snapshot = molecule->get_body(idx);
        auto cm = snapshot.get_cm();
        snapshot.translate(-cm);
        for (std::size_t j = 0; j < snapshot.size_symmetry(); ++j) {
            rigidbody->conformation->absolute_parameters.parameters[idx].symmetry_pars.emplace_back(snapshot.symmetry().get(j)->clone());
        }
        rigidbody->conformation->initial_conformation[idx] = std::move(snapshot);
        rigidbody->conformation->absolute_parameters.parameters[idx].translation = cm;
    }

    // the body count changed, so the histogram manager's per-body-indexed caches (sized for the old body count)
    // must be rebuilt from scratch, matching ConvertToSymmetryElement; the grid must be fully rebuilt too
    molecule->set_histogram_manager(std::make_unique<hist::PartialSymmetryManagerMT<true, false>>(molecule));
    rigidbody->molecule.clear_grid();
    rigidbody->refresh_grid();
}

std::vector<std::string> SplitElement::_valid_arguments() {
    return {};
}

std::unique_ptr<GenericElement> SplitElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("split", "Unexpected named argument \"" + args.named.begin()->first + "\".");}
    if (args.inlined.size() < 2) {
        throw except::parse_error("split", "Expected a body name followed by one or more residue sequence ids to split at.");
    }

    std::string body_name = args.inlined[0];
    std::vector<int> splits;
    splits.reserve(args.inlined.size()-1);
    for (std::size_t i = 1; i < args.inlined.size(); ++i) {
        const auto& tok = args.inlined[i];
        if (!utility::isinteger(tok)) {throw except::parse_error("split", "Invalid split index \"" + tok + "\"; expected an integer residue sequence id.");}
        splits.push_back(std::stoi(tok));
    }
    return std::make_unique<SplitElement>(owner->_get_sequencer(), body_name, std::move(splits));
}
