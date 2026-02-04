// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/OptimizeStepElement.h>
#include <rigidbody/sequencer/elements/OnImprovementElement.h>
#include <rigidbody/detail/Configuration.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Logging.h>

#include <iostream>

namespace ausaxs::rigidbody::sequencer {
    OptimizeStepElement::OptimizeStepElement(observer_ptr<LoopElement> owner) : LoopElement(owner, 1) {}

    OptimizeStepElement::~OptimizeStepElement() = default;

    void OptimizeStepElement::run() {
        logging::log("OptimizeStepElement::run: running optimization step " + std::to_string(LoopElement::global_counter));
        step_accepted = _get_sequencer()->_get_controller()->prepare_step();

        // Run all child elements after prepare_step() but before finish_step()
        // This allows them to see the attempted transformation
        for (auto& e : elements) {
            e->run();
        }

        _get_sequencer()->_get_controller()->finish_step();
        if (step_accepted) {
            logging::log("OptimizeStepElement::run: optimization step " + std::to_string(LoopElement::global_counter) + " accepted");
            if (settings::general::verbose) {
                std::cout << "Iteration " << LoopElement::global_counter << " of " << LoopElement::total_loop_count << std::endl;
                console::print_success("\tRigidBody::optimize: Accepted changes. New best chi2: " + std::to_string(_get_best_conf()->chi2));
            }
        }
    }
}