// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/Rigidbody.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>
#include <rigidbody/sequencer/elements/setup/SymmetryElement.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <settings/HistogramSettings.h>
#include <data/Molecule.h>
#include <data/Body.h>

#include <cassert>

using namespace ausaxs::rigidbody::sequencer;

namespace {
    using namespace ausaxs;

    // Place the symmetric copies at a sane distance from the original by seeding the
    // translation offset. CompositeSymmetry has no span of its own, so descend into its parts.
    void seed_offset(symmetry::ISymmetry& sym, double offset) {
        if (auto* comp = dynamic_cast<symmetry::CompositeSymmetry*>(&sym)) {
            seed_offset(*comp->inner, offset);
            seed_offset(*comp->outer, offset);
            return;
        }
        if (auto t = sym.span_translation(); !t.empty()) {*t.begin() = offset;}
    }

    // Mark the body's symmetry storage as optimizable. All supported symmetry types expose
    // both an offset and a frame orientation to the optimiser.
    void enable_optimization(symmetry::SymmetryStorage* storage) {
        auto* opt = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(storage);
        assert(opt != nullptr && "SymmetryElement: body symmetry storage is not optimizable.");
        opt->optimize_translate = true;
        opt->optimize_rot_axis = true;
    }
}

SymmetryElement::SymmetryElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<symmetry::type>& symmetry)
    : owner(owner)
{
    assert(names.size() == symmetry.size() && "SymmetryElement::SymmetryElement: The number of names and symmetries must be equal.");
    std::vector<std::unique_ptr<symmetry::ISymmetry>> symmetries;
    symmetries.reserve(symmetry.size());
    for (auto t : symmetry) {symmetries.emplace_back(symmetry::get(t));}
    _add(names, std::move(symmetries));
}

SymmetryElement::SymmetryElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<std::string>& symmetry_names)
    : owner(owner)
{
    assert(names.size() == symmetry_names.size() && "SymmetryElement::SymmetryElement: The number of names and symmetries must be equal.");
    std::vector<std::unique_ptr<symmetry::ISymmetry>> symmetries;
    symmetries.reserve(symmetry_names.size());
    for (const auto& name : symmetry_names) {symmetries.emplace_back(symmetry::create(name));}
    _add(names, std::move(symmetries));
}

void SymmetryElement::_add(const std::vector<std::string>& names, std::vector<std::unique_ptr<symmetry::ISymmetry>> symmetries) {
    assert(names.size() == symmetries.size() && "SymmetryElement::_add: The number of names and symmetries must be equal.");
    auto molecule = owner->_get_molecule();
    auto rigidbody = owner->_get_rigidbody();
    auto& setup = owner->setup();

    molecule->set_histogram_manager(std::make_unique<hist::PartialSymmetryManagerMT<true, false>>(molecule));
    for (unsigned int i = 0; i < names.size(); ++i) {
        auto index = setup._get_body_index(names[i]);
        int ibody = index.body;
        assert(index.replica == 0 && index.symmetry == -1 && "SymmetryElement::_add: The body name must refer to a base body (symmetry -1, replica 0).");

        // install the symmetry on the live body and the stored initial conformation; both
        // storages must be optimizable for the parameter optimiser to drive them
        enable_optimization(molecule->get_body(ibody).symmetry().get_obj());
        enable_optimization(rigidbody->conformation->initial_conformation[ibody].symmetry().get_obj());
        molecule->get_body(ibody).symmetry().add(symmetries[i]->clone());
        rigidbody->conformation->initial_conformation[ibody].symmetry().add(std::move(symmetries[i]));

        // add names for the symmetric bodies
        auto& name_map = setup._get_body_names();
        int isymmetry = molecule->get_body(ibody).size_symmetry()-1;
        assert(0 <= isymmetry && "SymmetryElement::_add: Inconsistent data structures.");
        if (int reps = molecule->get_body(ibody).symmetry().get(isymmetry)->repetitions(); reps == 1) { // single replica only: b1s1
            name_map.emplace(names[i] + "s" + std::to_string(isymmetry+1), detail::to_index(ibody, isymmetry, 1));
        } else { // multiple replicas, so include replica index in name: b1s1r1, b1s1r2, ...
            for (int j = 0; j < reps; ++j) {
                name_map.emplace(names[i] + "s" + std::to_string(isymmetry+1) + "r" + std::to_string(j+1), detail::to_index(ibody, isymmetry, j+1));
            }
        }

        // place the symmetry body at a sane distance from the original
        seed_offset(*molecule->get_body(ibody).symmetry().get(isymmetry), 2*molecule->get_Rg(false));
        seed_offset(*rigidbody->conformation->initial_conformation[ibody].symmetry().get(isymmetry), 2*molecule->get_Rg(false));
        rigidbody->conformation->absolute_parameters.parameters[ibody].symmetry_pars.emplace_back(
            molecule->get_body(ibody).symmetry().get(isymmetry)->clone()
        );
    }
    // Adding symmetry changes the body atom count, so the grid must be fully rebuilt
    rigidbody->molecule.clear_grid();
    rigidbody->refresh_grid();
}

SymmetryElement::~SymmetryElement() = default;

void SymmetryElement::run() {}

std::vector<std::string> SymmetryElement::_valid_arguments() {
    return {};
}

std::unique_ptr<GenericElement> SymmetryElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    auto rigidbody = owner->_get_rigidbody();

    // inline usage patterns
    if (!args.inlined.empty()) {
        if (!args.named.empty()) {throw except::parse_error("symmetry", "Cannot combine named and inline arguments.");}

        // usage pattern: [symmetry], allowed for single-body systems (e.g. "symmetry p2")
        if (args.inlined.size() == 1) {
            if (rigidbody->molecule.size_body() != 1) {throw except::parse_error("symmetry", "Could not determine which body to apply the symmetry to.");}
            return std::make_unique<SymmetryElement>(owner->_get_sequencer(), std::vector<std::string>{"b1"}, std::vector<std::string>{args.inlined[0]});
        }

        // usage pattern: [body] [symmetry] (e.g. "symmetry b1 p2")
        else if (args.inlined.size() == 2) {
            return std::make_unique<SymmetryElement>(owner->_get_sequencer(), std::vector<std::string>{args.inlined[0]}, std::vector<std::string>{args.inlined[1]});
        } else {
            throw except::parse_error("symmetry", "Too many inline arguments.");
        }
    }

    // named usage patterns; format is: 
    // symmetry {
    //    [body] [symmetry]
    //    [body] [symmetry]
    //    ...
    // }
    std::vector<std::string> symmetries;
    std::vector<std::string> names;
    for (auto& [name, value] : args.named) {
        if (value.size() != 1) {throw except::parse_error("symmetry", "Only one symmetry can be specified per body in each \"symmetry\" element.");}
        names.emplace_back(name);
        symmetries.emplace_back(value[0]);
    }
    return std::make_unique<SymmetryElement>(owner->_get_sequencer(), names, symmetries);
}