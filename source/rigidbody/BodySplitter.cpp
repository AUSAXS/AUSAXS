// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/BodySplitter.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/atoms/AtomMetadata.h>
#include <settings/MoleculeSettings.h>
#include <utility/Exceptions.h>
#include <io/pdb/PDBStructure.h>
#include <io/pdb/PDBAtom.h>
#include <io/pdb/PDBWater.h>
#include <io/Reader.h>

#include <algorithm>
#include <cassert>
#include <limits>
#include <optional>
#include <unordered_set>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::data;
using namespace ausaxs::io::pdb;

namespace {
    // Splitting by residue requires the residue ids to be retained as metadata. The operation is meaningless without
    // them, so the file-based overloads force retention on for the duration of the load rather than failing on a
    // configuration the caller has no reason to think about.
    struct residue_seq_guard {
        residue_seq_guard() : previous(settings::molecule::store_residue_seq) {settings::molecule::store_residue_seq = true;}
        ~residue_seq_guard() {settings::molecule::store_residue_seq = previous;}
        bool previous;
    };
}

std::vector<Body> BodySplitter::split(const Body& body, const std::vector<int>& splits) {
    const auto& metadata = body.get_metadata();
    if (!metadata || !metadata->residue_seq) {
        throw except::parse_error(
            "BodySplitter::split: Body has no residue sequence metadata. Enable settings::molecule::store_residue_seq "
            "before loading to make splitting by residue possible."
        );
    }
    const auto& atoms = body.get_atoms();
    const auto& resseq = *metadata->residue_seq;
    assert(resseq.size() == atoms.size() && "BodySplitter::split: residue_seq metadata is not parallel-indexed to the atom vector.");

    if (splits.empty()) {throw except::parse_error("BodySplitter::split: Expected at least one split index.");}
    if (atoms.empty()) {throw except::parse_error("BodySplitter::split: Cannot split an empty body.");}

    // note that residue ids may legitimately be negative, so the ids are used as lookup keys rather than as indices
    int max_id = std::numeric_limits<int>::min(), min_id = std::numeric_limits<int>::max();
    for (int r : resseq) {
        min_id = std::min(min_id, r);
        max_id = std::max(max_id, r);
    }

    // every split id must mark the first residue of a body that ends up with at least one atom in it. An empty body
    // has no observable presence beyond shifting all subsequent body indices, and its centre of mass is NaN (see
    // Body::get_cm, which divides by a total mass of zero), so it is rejected rather than produced.
    for (int id : splits) {
        if (id < min_id) {throw except::parse_error(
            "BodySplitter::split: Split " + std::to_string(id) + " "
            "smaller than lowest residue sequence id (" + std::to_string(min_id) + ")."
        );}
        if (max_id < id) {throw except::parse_error(
            "BodySplitter::split: Split " + std::to_string(id) + " "
            "larger than highest residue sequence id (" + std::to_string(max_id) + ")."
        );}
        if (id == min_id) {throw except::parse_error(
            "BodySplitter::split: Split " + std::to_string(id) + " "
            "is the lowest residue sequence id, which would leave the first body empty."
        );}
    }

    // each split id fires on the first atom carrying it and is consumed in the process, so a residue id occurring
    // again later in the body does not trigger a second split
    std::unordered_set<int> pending(splits.begin(), splits.end());
    if (pending.size() != splits.size()) {
        throw except::parse_error("BodySplitter::split: Duplicate split indices are not allowed, as they cannot each produce a body.");
    }

    auto slice = [&atoms, &metadata] (std::size_t begin, std::size_t end) -> Body {
        // the empty water vector is deliberate rather than incidental: it gives the body an explicit (currently empty)
        // hydration rather than the EmptyHydration null object, so that the resulting bodies can still be hydrated.
        // Any waters on the source body are dropped, as there is no well-defined assignment of them to the slices.
        Body b(std::vector<AtomFF>(atoms.begin()+begin, atoms.begin()+end), std::vector<data::Water>{});

        AtomMetadata m;
        bool any = false;
        if (metadata->backbone)    {m.backbone    = std::vector<backbone_t>(metadata->backbone->begin()+begin, metadata->backbone->begin()+end);   any = true;}
        if (metadata->residue_seq) {m.residue_seq = std::vector<int>(metadata->residue_seq->begin()+begin, metadata->residue_seq->begin()+end);     any = true;}
        if (metadata->occupancy)   {m.occupancy   = std::vector<float>(metadata->occupancy->begin()+begin, metadata->occupancy->begin()+end);       any = true;}
        if (any) {b.set_metadata(std::move(m));}
        return b;
    };

    std::vector<Body> bodies(splits.size()+1);
    std::size_t index_body = 0; // current index in the bodies vector
    std::size_t begin = 0;      // index in the atoms vector where the current body starts
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (pending.erase(resseq[i]) == 0) {continue;}
        bodies[index_body++] = slice(begin, i);
        begin = i;
    }
    bodies[index_body] = slice(begin, atoms.size());

    // an id left in `pending` never matched an atom, so it produced no body and left an unfilled entry at the tail.
    // Since the range check above bounds every id by an id that does occur, the only way to get here is a gap in the
    // residue numbering - common in structures with unresolved loops - which is almost certainly a mistake.
    if (!pending.empty()) {
        throw except::parse_error(
            "BodySplitter::split: Split " + std::to_string(*pending.begin()) + " does not match any residue in the body. "
            "The residue sequence ids likely contain a gap at this position."
        );
    }
    assert(index_body+1 == bodies.size() && "BodySplitter::split: not every body was filled.");
    return bodies;
}

Molecule BodySplitter::split(const io::File& input, const std::vector<int>& splits) {
    residue_seq_guard guard;
    return Molecule(split(Body(input), splits));
}

data::Molecule BodySplitter::split(const io::File& input) {
    io::pdb::PDBStructure data = io::Reader::read(input);
    if (settings::molecule::implicit_hydrogens) {data.add_implicit_hydrogens();}
    std::vector<PDBAtom>& atoms = data.atoms;

    std::vector<Body> bodies;
    auto begin = atoms.begin();
    char current_id = atoms[0].chainID;
    for (unsigned int i = 0; i < atoms.size(); i++) {
        if (atoms[i].chainID != current_id) {
            std::vector<PDBAtom> a(begin, atoms.begin() + i);
            auto reduced = PDBStructure(a, {}).reduced_representation();
            bodies.push_back(Body(reduced.atoms, reduced.waters));
            if (reduced.metadata) {bodies.back().set_metadata(std::move(*reduced.metadata));}
            begin = atoms.begin() + i;
            current_id = atoms[i].chainID;
        }
    }
    auto reduced = PDBStructure(std::vector<PDBAtom>(begin, atoms.end()), {}).reduced_representation();
    bodies.push_back(Body(reduced.atoms, reduced.waters));
    if (reduced.metadata) {bodies.back().set_metadata(std::move(*reduced.metadata));}
    return Molecule(bodies);
}
