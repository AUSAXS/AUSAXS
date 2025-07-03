// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/OptimizeStepElement.h>
#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/SaveElement.h>
#include <rigidbody/detail/BestConf.h>
#include <rigidbody/RigidBody.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>

#include <iostream>

namespace ausaxs::rigidbody::sequencer {
    OptimizeStepElement::OptimizeStepElement(LoopElement* owner) : LoopElementCallback(owner) {}

    OptimizeStepElement::~OptimizeStepElement() = default;

    void OptimizeStepElement::run() {
        if (owner->_get_sequencer()->_optimize_step()) {
            if (settings::general::verbose) {
                std::cout << "Iteration " << LoopElement::global_counter << " of " << LoopElement::total_loop_count << std::endl;
                console::print_success("\tRigidBody::optimize: Accepted changes. New best chi2: " + std::to_string(owner->_get_best_conf()->chi2));
            }

            for (auto& e : elements) {
                e->run();
            }
        }
    }
}