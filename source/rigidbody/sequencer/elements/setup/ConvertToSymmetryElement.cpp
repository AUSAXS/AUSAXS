// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/ConvertToSymmetryElement.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>
#include <rigidbody/sequencer/detail/SymmetryFit.h>
#include <rigidbody/sequencer/detail/BodyIndexOps.h>
#include <rigidbody/sequencer/detail/BodyNameRegistry.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/selection/SymmetryTargets.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <settings/GeneralSettings.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/Logging.h>

#include <algorithm>
#include <cassert>
#include <iostream>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    // Mark the body's symmetry storage as optimizable so the parameter optimiser can drive it.
    void enable_optimization(symmetry::SymmetryStorage* storage) {
        auto* opt = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(storage);
        assert(opt != nullptr && "ConvertToSymmetryElement: body symmetry storage is not optimizable.");
        opt->optimize_translate = true;
        opt->optimize_rot_axis = true;
    }
}

ConvertToSymmetryElement::ConvertToSymmetryElement(observer_ptr<Sequencer> owner, std::vector<int> bodies, const std::string& symmetry_name, double tolerance)
    : owner(owner)
{
    _convert(bodies, symmetry_name, tolerance);
}

ConvertToSymmetryElement::~ConvertToSymmetryElement() = default;

void ConvertToSymmetryElement::run() {}

void ConvertToSymmetryElement::_convert(const std::vector<int>& bodies, const std::string& symmetry_name, double tolerance) {
    detail::require_mutable_structure(owner, "convert_to_symmetry");

    auto molecule = owner->_get_molecule();
    auto rigidbody = owner->_get_rigidbody();
    auto& setup = owner->setup();

    // resolve and validate the requested symmetry type
    auto base_sym = symmetry::create(symmetry_name);
    if (!symmetry::is_optimizable(*base_sym)) {
        throw except::parse_error("convert_to_symmetry", "Unsupported symmetry \"" + symmetry_name + "\"; only point, cyclic and polyhedral symmetries can be fitted.");
    }

    if (bodies.size() != base_sym->repetitions() + 1) {
        throw except::parse_error("convert_to_symmetry",
            "Symmetry \"" + symmetry_name + "\" needs exactly " + std::to_string(base_sym->repetitions() + 1)
            + " bodies, but " + std::to_string(bodies.size()) + " were given.");
    }
    for (int b : bodies) {
        if (b < 0 || static_cast<std::size_t>(b) >= molecule->size_body()) {
            throw except::parse_error("convert_to_symmetry", "Body index out of range.");
        }
    }

    int primary = bodies.front();

    // gather the world-space atom coordinates of every participating body; correspondence is exact because they are copies of the same molecule 
    // (copies[k][i] is the image of copies[0][i])
    std::vector<std::vector<Vector3<double>>> copies;
    copies.reserve(bodies.size());
    for (int b : bodies) {
        const auto& atoms = molecule->get_body(b).get_atoms();
        if (atoms.size() != molecule->get_body(primary).get_atoms().size()) {
            throw except::parse_error("convert_to_symmetry", "All participating bodies must be copies of the same molecule (atom counts differ).");
        }
        std::vector<Vector3<double>> coords;
        coords.reserve(atoms.size());
        for (const auto& a : atoms) {coords.push_back(a.coordinates());}
        copies.push_back(std::move(coords));
    }
    Vector3<double> reference_cm = molecule->get_body(primary).get_cm();

    // fit the symmetry parameters to the assembly; the bodies (PDB chains) may be in any order, so let the fitter search over orderings, 
    // accepting the first that comes within tolerance
    auto fit = detail::fit_symmetry_best_order(*base_sym, reference_cm, copies, tolerance);
    logging::log("ConvertToSymmetryElement: fitted " + symmetry_name + " with residual RMSD " + std::to_string(fit.rmsd) + " Å.");
    if (settings::general::verbose) {
        std::cout << "\tFitted " << symmetry_name << " to " << bodies.size() << " bodies (residual RMSD "
                  << fit.rmsd << " Å)." << std::endl;
    }

    // the user asserts the symmetry type; a large residual means the assembly does not actually obey it, so we refuse rather than proceed 
    // with meaningless parameters
    if (tolerance < fit.rmsd) {
        throw except::parse_error("convert_to_symmetry",
            "The assembly does not match a \"" + symmetry_name + "\" symmetry (residual RMSD "
            + std::to_string(fit.rmsd) + " Å exceeds the tolerance of " + std::to_string(tolerance)
            + " Å). Check the symmetry type and the body order, or raise the tolerance.");
    }

    // install the fitted symmetry on the primary body (live molecule and stored initial conformation); both storages must be optimizable for 
    // the parameter optimiser to drive them
    enable_optimization(molecule->get_body(primary).symmetry().get_obj());
    enable_optimization(rigidbody->conformation->initial_conformation[primary].symmetry().get_obj());
    molecule->get_body(primary).symmetry().add(fit.symmetry->clone());
    rigidbody->conformation->initial_conformation[primary].symmetry().add(fit.symmetry->clone());
    rigidbody->conformation->absolute_parameters.parameters[primary].symmetry_pars.emplace_back(fit.symmetry->clone());

    // capture the new symmetry slot/replica count now, before body removal shifts the primary index
    int isymmetry = static_cast<int>(molecule->get_body(primary).size_symmetry()) - 1;
    int reps = static_cast<int>(molecule->get_body(primary).symmetry().get(isymmetry)->repetitions());

    // drop the now-redundant copy bodies; erase_bodies also reindexes every surviving body name
    std::vector<int> to_remove(bodies.begin() + 1, bodies.end());
    int new_primary = primary - static_cast<int>(std::count_if(to_remove.begin(), to_remove.end(), [primary](int b) {return b < primary;}));
    detail::erase_bodies(owner, std::move(to_remove));

    // register names for the primary body's newly-added symmetry copies
    for (int j = 0; j < reps; ++j) {setup._body_name_registry().add_replica(new_primary, isymmetry, j + 1);}

    // rebuild the (now symmetry-aware) histogram manager for the reduced body set; this also rebinds the body signallers. The grid must be rebuilt 
    // since the atom count changed.
    molecule->reset_histogram_manager();
    rigidbody->molecule.clear_grid();
    rigidbody->refresh_grid();
    rigidbody->symmetry_targets->invalidate(); // the set of declared symmetries changed
    rigidbody->constraints->invalidate();      // the redundant copies are gone, so the per-body constraint map is keyed by stale indices
}

std::vector<std::string> ConvertToSymmetryElement::_valid_arguments() {
    return {"type", "bodies", "tolerance"};
}

std::unique_ptr<GenericElement> ConvertToSymmetryElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    auto rigidbody = owner->_get_rigidbody();
    auto sequencer = owner->_get_sequencer();

    // inline form: "convert_to_symmetry c4" uses every loaded body, in load order
    if (!args.inlined.empty()) {
        if (!args.named.empty()) {throw except::parse_error("convert_to_symmetry", "Cannot combine named and inline arguments.");}
        if (args.inlined.size() != 1) {throw except::parse_error("convert_to_symmetry", "The inline form takes exactly one symmetry name.");}
        std::vector<int> bodies(rigidbody->molecule.size_body());
        for (int i = 0; i < static_cast<int>(bodies.size()); ++i) {bodies[i] = i;}
        return std::make_unique<ConvertToSymmetryElement>(sequencer, std::move(bodies), std::string{args.inlined[0]});
    }

    // block form: explicit "type" + "bodies" list, with an optional "tolerance"
    auto type_it = args.named.find("type");
    auto bodies_it = args.named.find("bodies");
    auto tolerance_it = args.named.find("tolerance");
    if (type_it == args.named.end() || bodies_it == args.named.end()) {
        throw except::parse_error("convert_to_symmetry", "The block form requires both a \"type\" and a \"bodies\" entry.");
    }
    std::size_t expected_keys = 2 + (tolerance_it != args.named.end() ? 1 : 0);
    if (args.named.size() != expected_keys) {throw except::parse_error("convert_to_symmetry", "Unexpected arguments; only \"type\", \"bodies\" and \"tolerance\" are allowed.");}
    if (type_it->second.size() != 1) {throw except::parse_error("convert_to_symmetry", "\"type\" takes exactly one symmetry name.");}

    double tolerance = ConvertToSymmetryElement::default_tolerance;
    if (tolerance_it != args.named.end()) {
        if (tolerance_it->second.size() != 1) {throw except::parse_error("convert_to_symmetry", "\"tolerance\" takes exactly one value.");}
        try {tolerance = std::stod(std::string{tolerance_it->second[0].str});}
        catch (const std::exception&) {throw except::parse_error("convert_to_symmetry", "\"tolerance\" must be a number.");}
    }

    std::string symmetry_name = type_it->second[0];
    std::vector<int> bodies;
    for (std::size_t i = 0; i < bodies_it->second.size(); ++i) {
        auto index = sequencer->setup()._get_body_index(std::string{bodies_it->second[i].str});
        if (index.symmetry != -1 || index.replica != 0) {
            throw except::parse_error("convert_to_symmetry", "Body names must refer to base bodies.");
        }
        bodies.push_back(index.body);
    }
    if (bodies.size() < 2) {throw except::parse_error("convert_to_symmetry", "At least two bodies are required.");}
    return std::make_unique<ConvertToSymmetryElement>(sequencer, std::move(bodies), symmetry_name, tolerance);
}
