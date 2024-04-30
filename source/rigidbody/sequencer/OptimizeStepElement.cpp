/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/OptimizeStepElement.h>
#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/SaveElement.h>
#include <rigidbody/detail/BestConf.h>
#include <rigidbody/RigidBody.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>

#include <iostream>

namespace rigidbody::sequencer {
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

    OptimizeStepElement& OptimizeStepElement::save_on_improvement(const io::File& path) {
        elements.push_back(std::make_unique<SaveElement>(owner, path));
        return *this;
    }
}