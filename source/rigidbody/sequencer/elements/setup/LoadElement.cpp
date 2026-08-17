// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/ArgumentHelper.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/setup/LoadElement.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>
#include <rigidbody/constraints/AttractorConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/constraints/RepellerConstraint.h>
#include <rigidbody/continuation/ContinuationState.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <utility/StringUtils.h>
#include <settings/GeneralSettings.h>

#include <algorithm>

using namespace ausaxs::rigidbody::sequencer;

LoadElement::~LoadElement() = default;

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& paths, const std::vector<std::string>& body_names) : owner(owner) {
    if (auto loc = paths[0].find("%"); loc != std::string::npos) {
        resolved_paths = load_wildcarded(paths[0]);
        rigidbody = std::make_unique<Rigidbody>(data::Molecule(resolved_paths));
    } else {
        resolved_paths = paths;
        std::transform(paths.begin(), paths.end(), resolved_paths.begin(), [this] (const std::string& path) {return lookup_file(path).first;});
        rigidbody = std::make_unique<Rigidbody>(data::Molecule(resolved_paths));
    }

    if (!body_names.empty() && body_names.size() != rigidbody->molecule.size_body()) {
        throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");
    }
    for (unsigned int i = 0; i < rigidbody->molecule.size_body(); ++i) {
        owner->setup()._body_name_registry().add_body(i, body_names.empty() ? std::string{} : body_names[i]);
    }
    owner->setup()._set_active_body(rigidbody.get());

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->molecule.size_body() << " bodies from " << paths.size() << " files." << std::endl;
    }
}

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<int>& splits, const std::vector<std::string>& body_names) 
    : owner(owner) 
{
    if (auto loc = path.find("%"); loc != std::string::npos) {
        resolved_paths = load_wildcarded(path);
        rigidbody = std::make_unique<Rigidbody>(data::Molecule(resolved_paths));
    } else {
        resolved_paths = {lookup_file(path).first};
        rigidbody = std::make_unique<Rigidbody>(rigidbody::BodySplitter::split(resolved_paths[0], splits));
    }

    if (!body_names.empty() && body_names.size() != rigidbody->molecule.size_body()) {
        throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");
    }
    for (unsigned int i = 0; i < rigidbody->molecule.size_body(); ++i) {
        owner->setup()._body_name_registry().add_body(i, body_names.empty() ? std::string{} : body_names[i]);
    }
    owner->setup()._set_active_body(rigidbody.get());

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->molecule.size_body() << " bodies from \"" << path << "\"." << std::endl;
    }
}

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<std::string>& body_names) : owner(owner) {
    resolved_paths = {lookup_file(path).first};
    rigidbody = std::make_unique<Rigidbody>(rigidbody::BodySplitter::split(resolved_paths[0]));
    if (rigidbody->molecule.size_body() <= 1) {
        throw std::runtime_error("LoadElement::LoadElement: Could not split \"" + path + "\" by chain, as it contains only one.");
    }

    if (!body_names.empty() && body_names.size() != rigidbody->molecule.size_body()) {
        throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");
    }
    for (unsigned int i = 0; i < rigidbody->molecule.size_body(); ++i) {
        owner->setup()._body_name_registry().add_body(i, body_names.empty() ? std::string{} : body_names[i]);
    }
    owner->setup()._set_active_body(rigidbody.get());

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->molecule.size_body() << " bodies from \"" << path << "\"." << std::endl;
    }
}

LoadElement::LoadElement(observer_ptr<Sequencer> owner, from_continuation_t, const std::string& path) : owner(owner) {
    auto [resolved, exists] = lookup_file(path);
    if (!exists) {throw std::runtime_error("LoadElement::LoadElement: Could not find continuation state file \"" + path + "\".");}
    resolved_paths = {resolved};

    auto state = continuation::read_continuation_state(io::File(resolved));
    rigidbody = std::make_unique<Rigidbody>(std::move(state.molecule));

    auto& names = owner->setup()._body_name_registry();
    for (unsigned int i = 0; i < rigidbody->molecule.size_body(); ++i) {
        names.add_body(i, i < state.body_names.size() ? state.body_names[i] : std::string{});

        // replicas are registered here rather than by a symmetry element, since a continued run declares no symmetry
        // elements of its own — the symmetries arrive already attached to the restored bodies.
        const auto& symmetries = rigidbody->molecule.get_body(i).symmetry().get();
        for (int isymmetry = 0; isymmetry < static_cast<int>(symmetries.size()); ++isymmetry) {
            for (int replica = 1; replica <= static_cast<int>(symmetries[isymmetry]->repetitions()); ++replica) {
                names.add_replica(i, isymmetry, replica);
            }
        }
    }

    // constraints are rebuilt from their stored description rather than re-derived: see constraints::restore_t
    for (const auto& c : state.constraints) {
        std::unique_ptr<constraints::Constraint> constraint;
        switch (c.kind) {
            using Kind = continuation::ContinuationConstraint::Kind;
            case Kind::bond:
                constraint = std::make_unique<constraints::DistanceConstraintBond>(
                    constraints::restore, &rigidbody->molecule, c.ibody1, c.iatom1, c.ibody2, c.iatom2, c.isym1, c.isym2, c.d_target);
                break;
            case Kind::cm:
                constraint = std::make_unique<constraints::DistanceConstraintCM>(
                    constraints::restore, &rigidbody->molecule, c.ibody1, c.iatom1, c.ibody2, c.iatom2, c.isym1, c.isym2, c.d_target);
                break;
            case Kind::attractor:
                constraint = std::make_unique<constraints::AttractorConstraint>(
                    constraints::restore, &rigidbody->molecule, c.ibody1, c.iatom1, c.ibody2, c.iatom2, c.isym1, c.isym2, c.d_target);
                break;
            case Kind::repeller:
                constraint = std::make_unique<constraints::RepellerConstraint>(
                    constraints::restore, &rigidbody->molecule, c.ibody1, c.iatom1, c.ibody2, c.iatom2, c.isym1, c.isym2, c.d_target);
                break;
            case Kind::atom:
                constraint = std::make_unique<constraints::DistanceConstraintAtom>(
                    constraints::restore, &rigidbody->molecule, c.ibody1, c.iatom1, c.ibody2, c.iatom2, c.isym1, c.isym2, c.d_target);
                break;
        }
        rigidbody->constraints->add_constraint(std::move(constraint));
    }

    owner->setup()._set_active_body(rigidbody.get());

    if (settings::general::verbose) {
        std::cout << "\tResumed " << rigidbody->molecule.size_body() << " bodies and " << state.constraints.size()
                  << " constraints from \"" << resolved << "\"." << std::endl;
    }
}

/**
 * @brief Looks up a file relative to both the current working directory and the configuration file directory.
 * @return A pair containing the resolved file path and a boolean indicating whether the file exists.
 */
std::pair<std::string, bool> LoadElement::lookup_file(const std::string& path) {
    io::File file(path);
    if (file.exists()) {return {file, true};}

    auto config_folder = owner->setup()._get_config_folder();
    io::File relative(config_folder, file.stem(), file.extension());
    if (relative.exists()) {return {relative, true};}

    return {file, false};
}

std::vector<std::string> LoadElement::load_wildcarded(const std::string& path) {
    static auto zero_pad_string = [] (int val, unsigned int pad) -> std::string {
        std::string s = std::to_string(val);
        if (s.size() < pad) {
            s.insert(0, pad - s.size(), '0');
        }
        return s;
    };
    std::vector<std::string> wildcarded_files;

    auto loc = path.find("%");
    int start = loc, end = loc;
    while (path[++end] == '%') {if (100 < end - start) {throw std::runtime_error("LoadElement::LoadElement: The maximum number of consecutive '%' characters is 100.");}}
    int counter = 0;

    // file numbered zero may or may not exist
    io::File file = path.substr(0, start) + zero_pad_string(counter, end - start) + path.substr(end);
    if (auto res = lookup_file(file); res.second) { // check filename padded with zeros
        wildcarded_files.push_back(res.first);
    } else if (1 < end - start) {   // check filename without padding
        file = path.substr(0, start) + std::to_string(counter) + path.substr(end);
        if (res = lookup_file(file); res.second) {wildcarded_files.push_back(res.first);}
    }
    ++counter;

    // check for wildcarded files with padding
    file = path.substr(0, start) + zero_pad_string(counter, end - start) + path.substr(end);
    auto res = lookup_file(file);
    while (res.second) {
        wildcarded_files.push_back(res.first);
        file = path.substr(0, start) + zero_pad_string(++counter, end - start) + path.substr(end);
        res = lookup_file(file);
    }

    // check for wildcarded files without padding
    if (wildcarded_files.size() <= 1) { // only check if no files were found with padding
        file = path.substr(0, start) + std::to_string(counter) + path.substr(end);
        res = lookup_file(file);
        while (res.second) {
            wildcarded_files.push_back(file.path());
            file = path.substr(0, start) + std::to_string(++counter) + path.substr(end);
            res = lookup_file(file);
        }
    }

    if (wildcarded_files.empty()) {throw std::runtime_error("LoadElement::LoadElement: No files found matching the wildcarded path.");}
    return wildcarded_files;
}

void LoadElement::run() {
    owner->_get_sequencer()->_set_rigidbody(rigidbody.get());
}

namespace {
    enum class Args {paths, splits, names, saxs, continuation};
    std::unordered_map<Args, std::vector<std::string>> args_map = {
        {Args::paths, {"pdb"}},
        {Args::saxs, {"saxs"}},
        {Args::splits, {"split"}},
        {Args::names, {"names", "name"}},
        {Args::continuation, {"continue"}}
    };
}

std::vector<std::string> LoadElement::_valid_arguments() {
    static auto map = detail::get_arg_names<Args>(args_map);
    return map;
}

std::unique_ptr<GenericElement> LoadElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    enum class Args {paths, splits, names, saxs, continuation};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::paths, {"pdb"}},
        {Args::saxs, {"saxs"}},
        {Args::splits, {"split"}},
        {Args::names, {"names", "name"}},
        {Args::continuation, {"continue"}}
    };

    auto pdb = args.get<std::vector<std::string>>(valid_args[Args::paths]);
    auto saxs = args.get<std::string>(valid_args[Args::saxs]);
    auto names = args.get<std::vector<std::string>>(valid_args[Args::names]);
    auto split = args.get<std::vector<std::string>>(valid_args[Args::splits]);
    auto resume = args.get<std::string>(valid_args[Args::continuation]);

    if (!args.inlined.empty()) {throw except::parse_error("load", "Unexpected inline arguments.");}
    if (!saxs.found) {throw except::parse_error("load", "Missing required argument \"saxs\".");}

    if (resume.found) {
        // a continuation state already fixes the body decomposition and their names, so the arguments that would set
        // those up are rejected rather than silently ignored
        if (pdb.found)   {throw except::parse_error("load", "\"continue\" and \"pdb\" are mutually exclusive.");}
        if (split.found) {throw except::parse_error("load", "\"split\" cannot be combined with \"continue\"; the continuation state already carries the body decomposition.");}
        if (names.found) {throw except::parse_error("load", "\"names\" cannot be combined with \"continue\"; the continuation state already carries the body names.");}
        owner->_get_sequencer()->setup()._set_saxs_path(io::ExistingFile(saxs.value));
        return std::make_unique<LoadElement>(owner->_get_sequencer(), from_continuation, resume.value);
    }
    if (!pdb.found) {throw except::parse_error("load", "Missing required argument \"path\".");}

    owner->_get_sequencer()->setup()._set_saxs_path(io::ExistingFile(saxs.value));
    if (split.found) {
        if (split.value.size() == 1 && split.value[0] == "chain") {
            if (pdb.value.size() != 1) {throw except::parse_error("load", "Chain splitting can only be used with a single path.");}
            return std::make_unique<LoadElement>(owner->_get_sequencer(), pdb.value[0], names.value);
        } else {
            std::vector<int> splits;
            for (const auto& s : split.value) {
                if (!utility::isinteger(s)) {
                    throw except::parse_error("load", "Invalid argument for \"split\": \"" + s + "\". Expected \"chain\" or a list of positive integers.");
                }
                splits.push_back(std::stoi(s));
            }
            return std::make_unique<LoadElement>(owner->_get_sequencer(), pdb.value[0], splits, names.value);
        }
    }
    return std::make_unique<LoadElement>(owner->_get_sequencer(), pdb.value, names.value);
}