// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/ConstraintElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/constraints/ConstraintFactory.h>
#include <rigidbody/Rigidbody.h>

using namespace ausaxs::rigidbody::sequencer;

ConstraintElement::ConstraintElement(observer_ptr<Sequencer> owner, std::unique_ptr<rigidbody::constraints::Constraint> constraint) {
    owner->_get_rigidbody()->constraints->add_constraint(std::move(constraint));
}

void ConstraintElement::run() {}

namespace {
    enum class ConstraintChoice {
        SpecificAtoms,
        Bond,
        BodyCM,
        BodyCMAttractor,
        BodyCMRepeller,
        Unknown
    };
}

void ConstraintElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    enum class Args {body1, body2, iatom1, iatom2, type, distance};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::body1, {"first", "body1"}},
        {Args::body2, {"second", "body2"}},
        {Args::iatom1, {"iatom1", "iatom_1", "i1"}},
        {Args::iatom2, {"iatom2", "iatom_2", "i2"}},
        {Args::type, {"type", "kind"}},
        {Args::distance, {"distance", "dist", "d"}}
    };

    auto body1 = args.get<std::string>(valid_args[Args::body1]);
    auto body2 = args.get<std::string>(valid_args[Args::body2]);
    auto type  = args.get<std::string>(valid_args[Args::type]);
    auto iatom1 = args.get<int>(valid_args[Args::iatom1]);
    auto iatom2 = args.get<int>(valid_args[Args::iatom2]);
    auto distance = args.get<double>(valid_args[Args::distance]);

    static auto get_constraint_choice = [] (std::string_view line) {
        if (line == "bond") {return ConstraintChoice::Bond;}
        if (line == "distance") {return ConstraintChoice::SpecificAtoms;}
        if (line == "cm" || line == "center_of_mass") {return ConstraintChoice::BodyCM;}
        if (line == "attract") {return ConstraintChoice::BodyCMAttractor;}
        if (line == "repel") {return ConstraintChoice::BodyCMRepeller;}
        throw except::parse_error("constraint", "Unknown choice \"" + std::string(line) + "\"");
    };

    if (!args.inlined.empty()) {throw except::parse_error("constraint", "Unexpected inline arguments.");}
    if (!body1.found) {throw except::parse_error("constraint", "Missing required argument \"body1\".");}
    if (!body2.found) {throw except::parse_error("constraint", "Missing required argument \"body2\".");}
    if (!type.found) {throw except::parse_error("constraint", "Missing required argument \"type\".");}

    std::unique_ptr<constraints::Constraint> constraint;
    ConstraintChoice type_enum = get_constraint_choice(type.value);
    switch (type_enum) {
        case ConstraintChoice::Bond:
            if (iatom1.found || iatom2.found || distance.found) {
                throw except::parse_error("constraint", "Constraint of type \"bond\" cannot be provided with arguments \"iatom1\", \"iatom2\" or \"distance\".");
            }
            constraint = factory::create_constraint_bond(
                owner->_get_molecule(),
                owner->_get_sequencer()->setup()._get_body_index(body1.value),
                owner->_get_sequencer()->setup()._get_body_index(body2.value)
            );
            break;

        case ConstraintChoice::BodyCM:
            if (iatom1.found || iatom2.found || distance.found) {
                throw except::parse_error("constraint", "Constraint of type \"cm\" cannot be provided with arguments \"iatom1\", \"iatom2\" or \"distance\".");
            }
            constraint = factory::create_constraint_cm(
                owner->_get_molecule(),
                owner->_get_sequencer()->setup()._get_body_index(body1.value),
                owner->_get_sequencer()->setup()._get_body_index(body2.value)
            );
            break;

        case ConstraintChoice::BodyCMAttractor:
            if (iatom1.found || iatom2.found) {
                throw except::parse_error("constraint", "Constraint of type \"attract\" cannot be provided with arguments \"iatom1\" or \"iatom2\".");
            }
            if (!distance.found) {throw except::parse_error("constraint", "Constraint of type \"attract\" requires argument \"distance\".");}
            constraint = factory::create_constraint_attractor(
                owner->_get_molecule(),
                owner->_get_sequencer()->setup()._get_body_index(body1.value),
                owner->_get_sequencer()->setup()._get_body_index(body2.value),
                distance.value
            );
            break;

        case ConstraintChoice::BodyCMRepeller:
            if (iatom1.found || iatom2.found) {
                throw except::parse_error("constraint", "Constraint of type \"repel\" cannot be provided with arguments \"iatom1\" or \"iatom2\".");
            }
            if (!distance.found) {throw except::parse_error("constraint", "Constraint of type \"repel\" requires argument \"distance\".");}
            constraint = factory::create_constraint_repeller(
                owner->_get_molecule(),
                owner->_get_sequencer()->setup()._get_body_index(body1.value),
                owner->_get_sequencer()->setup()._get_body_index(body2.value),
                distance.value
            );
            break;

        case ConstraintChoice::SpecificAtoms:
            if (distance.found) {throw except::parse_error("constraint", "Constraint of type \"specific_atoms\" cannot be provided with argument \"distance\".");}
            if (!(iatom1.found && iatom2.found)) {throw except::parse_error("constraint", "Constraint of type \"specific_atoms\" requires arguments \"iatom1\" and \"iatom2\".");}
            constraint = factory::create_constraint(
                owner->_get_molecule(),
                owner->_get_sequencer()->setup()._get_body_index(body1.value),
                owner->_get_sequencer()->setup()._get_body_index(body2.value),
                iatom1.value,
                iatom2.value
            );
            break;

        default:
            throw except::parse_error("constraint", "Unknown choice \"" + type.value + "\".");
    }
    owner->_get_rigidbody()->constraints->add_constraint(std::move(constraint));
}