// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/LoadElementWithMetadata.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/SetupElement.h>
#include <rigidbody/sequencer/detail/ArgumentHelper.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/Reader.h>
#include <io/File.h>
#include <io/pdb/PDBStructure.h>
#include <io/pdb/PDBAtom.h>
#include <constants/ConstantsFwd.h>
#include <utility/StringUtils.h>

#include <algorithm>
#include <stdexcept>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

LoadElementWithMetadata::LoadElementWithMetadata(observer_ptr<Sequencer> owner, const std::vector<std::string>& paths, const std::vector<std::string>& body_names)
    : LoadElement(owner, paths, body_names) {capture_metadata();}

LoadElementWithMetadata::LoadElementWithMetadata(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<int>& split, const std::vector<std::string>& body_names)
    : LoadElement(owner, path, split, body_names) {capture_metadata();}

LoadElementWithMetadata::LoadElementWithMetadata(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<std::string>& body_names)
    : LoadElement(owner, path, body_names) {capture_metadata();}

LoadElementWithMetadata::~LoadElementWithMetadata() = default;

void LoadElementWithMetadata::capture_metadata() {
    residue_seq.clear();
    is_ca.clear();

    // re-read the same files; loading discards hydrogens and stores hydration separately, so
    // these protein atoms are in 1:1 correspondence (same count and order) with the molecule
    for (const auto& path : resolved_paths) {
        auto structure = io::Reader::read(io::File(path));
        residue_seq.reserve(residue_seq.size() + structure.atoms.size());
        is_ca.reserve(is_ca.size() + structure.atoms.size());
        for (const auto& atom : structure.atoms) {
            std::string name = atom.name;
            name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
            residue_seq.push_back(atom.resSeq);
            // a Cα is an atom named "CA" of element carbon (excludes calcium ions, also named "CA")
            is_ca.push_back(name == "CA" && atom.element == constants::atom_t::C ? 1 : 0);
        }
    }

    std::size_t n_molecule = 0;
    for (const auto& body : rigidbody->molecule.get_bodies()) {n_molecule += body.size_atom();}
    if (n_molecule != residue_seq.size()) {
        throw std::runtime_error(
            "LoadElementWithMetadata: atom-count mismatch between the molecule (" + std::to_string(n_molecule) +
            ") and the re-read PDB metadata (" + std::to_string(residue_seq.size()) +
            "); cannot align the structure preview."
        );
    }
}

std::unique_ptr<GenericElement> LoadElementWithMetadata::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    // mirrors LoadElement::_parse, but constructs the metadata-capturing variant instead
    enum class Args {paths, splits, names, saxs};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::paths, {"pdb"}},
        {Args::saxs, {"saxs"}},
        {Args::splits, {"split"}},
        {Args::names, {"names", "name"}}
    };

    auto pdb = args.get<std::vector<std::string>>(valid_args[Args::paths]);
    auto saxs = args.get<std::string>(valid_args[Args::saxs]);
    auto names = args.get<std::vector<std::string>>(valid_args[Args::names]);
    auto split = args.get<std::vector<std::string>>(valid_args[Args::splits]);

    if (!args.inlined.empty()) {throw except::parse_error("load", "Unexpected inline arguments.");}
    if (!pdb.found) {throw except::parse_error("load", "Missing required argument \"path\".");}
    if (!saxs.found) {throw except::parse_error("load", "Missing required argument \"saxs\".");}

    owner->_get_sequencer()->setup()._set_saxs_path(io::ExistingFile(saxs.value));
    if (split.found) {
        // strip trailing commas from each token to allow "split 15, 106, 206" style
        for (auto& s : split.value) {
            if (!s.empty() && s.back() == ',') { s.pop_back(); }
        }
        if (split.value.size() == 1 && split.value[0] == "chain") {
            if (pdb.value.size() != 1) {throw except::parse_error("load", "Chain splitting can only be used with a single path.");}
            return std::make_unique<LoadElementWithMetadata>(owner->_get_sequencer(), pdb.value[0], names.value);
        } else {
            std::vector<int> splits;
            for (const auto& s : split.value) {
                if (!utility::isinteger(s)) {
                    throw except::parse_error("load", "Invalid argument for \"split\": \"" + s + "\". Expected \"chain\" or a list of positive integers.");
                }
                splits.push_back(std::stoi(s));
            }
            return std::make_unique<LoadElementWithMetadata>(owner->_get_sequencer(), pdb.value[0], splits, names.value);
        }
    }
    return std::make_unique<LoadElementWithMetadata>(owner->_get_sequencer(), pdb.value, names.value);
}
