// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/OptimizeStepElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/elements/SaveElement.h>
#include <rigidbody/detail/Configuration.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Logging.h>

#include <iostream>

namespace ausaxs::rigidbody::sequencer {
    OptimizeStepElement::OptimizeStepElement(LoopElement* owner) : LoopElementCallback(owner) {}

    OptimizeStepElement::~OptimizeStepElement() = default;

    void OptimizeStepElement::run() {
        logging::log("OptimizeStepElement::run: running optimization step " + std::to_string(LoopElement::global_counter));
        owner->_get_sequencer()->_get_controller()->prepare_step();

        for (auto& e : elements) {
            e->run();
        }

        if (owner->_get_sequencer()->_get_controller()->finish_step()) {
            logging::log("OptimizeStepElement::run: optimization step " + std::to_string(LoopElement::global_counter) + " accepted");
            if (settings::general::verbose) {
                std::cout << "Iteration " << LoopElement::global_counter << " of " << LoopElement::total_loop_count << std::endl;
                console::print_success("\tRigidBody::optimize: Accepted changes. New best chi2: " + std::to_string(owner->_get_best_conf()->chi2));
            }
        }
    }
}