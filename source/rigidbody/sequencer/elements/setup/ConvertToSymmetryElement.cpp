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
#include <map>

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

    // A contiguous run of atoms sharing a residue id, keyed by that id together with its occurrence count so that a body holding the same id more than
    // once - a merge of several chains, each numbered from its own start - still matches its copies.
    struct ResidueRun {
        std::pair<int, int> key;
        std::size_t begin, size;
    };

    // Split a body into its residue runs, in atom order. Empty if the body carries no residue metadata to split on.
    std::vector<ResidueRun> residue_runs(const data::Body& body) {
        const auto& metadata = body.get_metadata();
        if (!metadata || !metadata->residue_seq) {return {};}
        const auto& seq = *metadata->residue_seq;
        assert(seq.size() == body.size_atom() && "residue_runs: metadata is not parallel-indexed to the atom vector.");

        std::vector<ResidueRun> runs;
        std::map<int, int> occurrences;
        for (std::size_t i = 0; i < seq.size();) {
            std::size_t j = i;
            while (j < seq.size() && seq[j] == seq[i]) {++j;}
            runs.push_back({{seq[i], occurrences[seq[i]]++}, i, j - i});
            i = j;
        }
        return runs;
    }

    // Reduce a body to its first n atoms, keeping the optional metadata parallel-indexed to the shortened atom vector.
    void truncate(data::Body& body, std::size_t n) {
        assert(n <= body.size_atom() && "ConvertToSymmetryElement::truncate: cannot grow a body.");
        body.get_atoms().resize(n);
        if (!body.get_metadata()) {return;}
        data::AtomMetadata metadata = *body.get_metadata();
        if (metadata.backbone)    {metadata.backbone->resize(n);}
        if (metadata.residue_seq) {metadata.residue_seq->resize(n);}
        if (metadata.occupancy)   {metadata.occupancy->resize(n);}
        body.set_metadata(std::move(metadata));
    }
}

ConvertToSymmetryElement::ConvertToSymmetryElement(observer_ptr<Sequencer> owner, std::vector<int> bodies, const std::string& symmetry_name, double tolerance)
    : owner(owner)
{
    _convert(bodies, symmetry_name, tolerance);
}

ConvertToSymmetryElement::~ConvertToSymmetryElement() = default;

void ConvertToSymmetryElement::run() {}

std::vector<std::vector<Vector3<double>>> ConvertToSymmetryElement::_split_into_copies(int primary, std::size_t copies_wanted, const std::string& symmetry_name) {
    auto molecule = owner->_get_molecule();
    auto rigidbody = owner->_get_rigidbody();
    auto& body = molecule->get_body(primary);

    // The copies of such a structure are laid out sequentially in the source file (each is a chain, or a contiguous run of residues within one), so equal
    // contiguous chunks of the atom vector recover them. Unequal chunks mean the copies are not identical, which the fit requires; a wrongly-guessed
    // decomposition instead surfaces as a large fit residual and is rejected by the tolerance check in _convert.
    std::size_t total = body.size_atom();
    if (total == 0 || total % copies_wanted != 0) {
        throw except::parse_error("convert_to_symmetry",
            "A single body was given, so it must itself be the assembly: \"" + symmetry_name + "\" requires it to split into "
            + std::to_string(copies_wanted) + " equally-sized copies, but its " + std::to_string(total)
            + " atoms do not divide evenly among them. Either the structure is not a whole number of copies, or its copies are not identical.");
    }
    std::size_t n = total/copies_wanted;

    std::vector<std::vector<Vector3<double>>> copies;
    copies.reserve(copies_wanted);
    for (std::size_t k = 0; k < copies_wanted; ++k) {
        std::vector<Vector3<double>> coords;
        coords.reserve(n);
        for (std::size_t i = k*n; i < (k+1)*n; ++i) {coords.push_back(body.get_atom(i).coordinates());}
        copies.push_back(std::move(coords));
    }

    // reduce the body to the leading copy, which the fitted symmetry will regenerate the rest from. The stored initial conformation is parallel-indexed to
    // the live body, so it is truncated identically; re-centring it on the origin preserves the invariant relied on when a transform rebuilds the body from
    // it (see TransformStrategy::apply), with the body's cm - which is where the new position lands - restored as its absolute translation.
    truncate(body, n);
    auto& initial = rigidbody->conformation->initial_conformation[primary];
    truncate(initial, n);
    initial.translate(-initial.get_cm());
    rigidbody->conformation->absolute_parameters.parameters[primary].translation = body.get_cm();

    logging::log("ConvertToSymmetryElement: split the single body into " + std::to_string(copies_wanted) + " copies of " + std::to_string(n) + " atoms each.");
    if (settings::general::verbose) {
        std::cout << "\tSplit the single body into " << copies_wanted << " copies of " << n << " atoms each." << std::endl;
    }
    return copies;
}

std::vector<std::vector<Vector3<double>>> ConvertToSymmetryElement::_gather_copies(const std::vector<int>& bodies) {
    auto molecule = owner->_get_molecule();
    int primary = bodies.front();
    auto body_name = [this](int b) {return owner->setup()._body_name_registry().base_body_names().at(b);};

    // Copies of the same molecule are not necessarily modelled to the same extent: a disordered terminus or loop is routinely resolved in some chains and not
    // in others, leaving their atom vectors differing in both length and content. Since the fit needs a correspondence rather than every atom, it runs on the
    // residues every copy has in common. Matching on residue identity - rather than trimming to the shortest atom vector - is what keeps the correspondence
    // right when the unmodelled stretch is interior instead of terminal, where every atom past the gap would otherwise pair up with the wrong one.
    std::vector<std::vector<ResidueRun>> runs;
    runs.reserve(bodies.size());
    for (int b : bodies) {runs.push_back(residue_runs(molecule->get_body(b)));}

    for (std::size_t k = 0; k < bodies.size(); ++k) {
        if (runs[k].empty()) {
            throw except::parse_error("convert_to_symmetry",
                "Body \"" + body_name(bodies[k]) + "\" carries no residue metadata, which the copies are matched up by. This should not be reachable: "
                "rigidbody refinement retains it unconditionally (see settings::molecule::store_residue_seq)."
            );
        }
    }

    std::vector<std::vector<Vector3<double>>> copies;
    copies.reserve(bodies.size());

    // Tally every residue over the participating bodies, keeping those all of them hold with the same number of atoms. A residue modelled to differing extents
    // (a partial side chain) offers no usable correspondence either, so it is dropped along with the ones that are missing outright.
    struct Tally {int bodies; std::size_t size;};
    std::map<std::pair<int, int>, Tally> tally;
    for (const auto& body_runs : runs) {
        for (const auto& run : body_runs) {
            auto [it, inserted] = tally.try_emplace(run.key, Tally{0, run.size});
            if (it->second.size != run.size) {it->second.bodies = -1;} // sticky: an atom count that differs anywhere disqualifies the residue everywhere
            else if (0 <= it->second.bodies) {++it->second.bodies;}
        }
    }

    for (std::size_t k = 0; k < bodies.size(); ++k) {
        const auto& atoms = molecule->get_body(bodies[k]).get_atoms();
        std::vector<Vector3<double>> coords;
        for (const auto& run : runs[k]) {
            if (tally.at(run.key).bodies != static_cast<int>(bodies.size())) {continue;}
            for (std::size_t i = run.begin; i < run.begin + run.size; ++i) {coords.push_back(atoms[i].coordinates());}
        }
        copies.push_back(std::move(coords));
    }
    assert(
        std::all_of(copies.begin(), copies.end(), [&copies] (const auto& c) {return c.size() == copies[0].size();})
        && "_gather_copies: the shared residues must contribute equally to every copy."
    );

    if (copies[0].empty()) {
        throw except::parse_error("convert_to_symmetry",
            "The participating bodies share no residues, so there is nothing to fit the symmetry to. They are either not copies of the same molecule, or "
            "their residue numbering does not agree."
        );
    }

    // the fit runs on the shared subset, but the primary body is kept whole and is what the fitted symmetry goes on to replicate
    std::size_t dropped = molecule->get_body(primary).size_atom() - copies[0].size();
    if (0 < dropped) {
        logging::log("ConvertToSymmetryElement: the copies are modelled to differing extents; fitted on the "
            + std::to_string(copies[0].size()) + " atoms they share, ignoring " + std::to_string(dropped) + " of the primary body's.");
        if (settings::general::verbose) {
            std::cout << "\tThe copies are modelled to differing extents; fitting on the " << copies[0].size()
                      << " atoms they share (ignoring " << dropped << " of the primary body's)." << std::endl;
        }
    }
    return copies;
}

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

    std::size_t expected = base_sym->repetitions() + 1;
    // a single body cannot be a set of copies, so it must be the assembled structure itself; the decomposition is then ours to find rather than the user's
    // to supply (see the auto-split branch below)
    bool auto_split = bodies.size() == 1 && 1 < expected;
    if (!auto_split && bodies.size() != expected) {
        throw except::parse_error("convert_to_symmetry",
            "Symmetry \"" + symmetry_name + "\" needs exactly " + std::to_string(expected)
            + " bodies, but " + std::to_string(bodies.size()) + " were given."
        );
    }
    for (int b : bodies) {
        if (b < 0 || static_cast<std::size_t>(b) >= molecule->size_body()) {
            throw except::parse_error("convert_to_symmetry", "Body index out of range.");
        }
    }

    int primary = bodies.front();

    // gather the world-space atom coordinates of every participating copy, index-parallel so that copies[k][i] is the image of copies[0][i]
    std::vector<std::vector<Vector3<double>>> copies;
    copies.reserve(expected);
    if (auto_split) {
        copies = _split_into_copies(primary, expected, symmetry_name);
    } else {
        copies = _gather_copies(bodies);
    }
    // in the auto-split case the body has already been reduced to the reference copy, so this is the reference's own cm either way. It is the whole body's
    // cm even when the fit only saw a subset of its atoms: the cm is the anchor the fitted symmetry is expressed relative to. 
    // The subset only determines the transform between copies, which is independent of the anchor.
    Vector3<double> reference_cm = molecule->get_body(primary).get_cm();

    // fit the symmetry parameters to the assembly; the copies (PDB chains, or the chunks of an auto-split body) may be in any order, so let the fitter
    // search over orderings, accepting the first that comes within tolerance
    auto fit = detail::fit_symmetry_best_order(*base_sym, reference_cm, copies, tolerance);
    logging::log("ConvertToSymmetryElement: fitted " + symmetry_name + " with residual RMSD " + std::to_string(fit.rmsd) + " Å.");
    if (settings::general::verbose) {
        std::cout << "\tFitted " << symmetry_name << " to " << copies.size() << " copies (residual RMSD "
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

    // drop the now-redundant copy bodies; erase_bodies also reindexes every surviving body name. An auto-split has no such bodies to drop - the copies were
    // only ever slices of the primary - so it keeps both its index and its name.
    std::vector<int> to_remove(bodies.begin() + 1, bodies.end());
    int new_primary = primary - static_cast<int>(std::count_if(to_remove.begin(), to_remove.end(), [primary](int b) {return b < primary;}));
    if (!to_remove.empty()) {detail::erase_bodies(owner, std::move(to_remove));}

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
    if (bodies.empty()) {throw except::parse_error("convert_to_symmetry", "At least one body is required.");}
    return std::make_unique<ConvertToSymmetryElement>(sequencer, std::move(bodies), symmetry_name, tolerance);
}
