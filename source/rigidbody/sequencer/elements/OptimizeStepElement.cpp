// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/OptimizeStepElement.h>
#include <rigidbody/sequencer/elements/OnImprovementElement.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Logging.h>

using namespace ausaxs::rigidbody::sequencer;
OptimizeStepElement::OptimizeStepElement(observer_ptr<LoopElement> owner) : LoopElement(owner, 1) {}

OptimizeStepElement::~OptimizeStepElement() = default;

void OptimizeStepElement::run() {
    step_accepted = _get_sequencer()->_get_controller()->prepare_step();

    // Run all child elements after prepare_step() but before finish_step()
    // This allows them to see the attempted transformation
    for (auto& e : elements) {
        e->run();
    }

    _get_sequencer()->_get_controller()->finish_step();
}

std::unique_ptr<GenericElement> OptimizeStepElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("optimize_step", "Unexpected named argument.");}
    if (!args.inlined.empty()) {throw except::parse_error("optimize_step", "Unexpected inline argument.");}
    return std::make_unique<OptimizeStepElement>(owner);
}