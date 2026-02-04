// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/LoopElement.h>

namespace ausaxs::rigidbody::sequencer {
    class OnImprovementElement;

    /**
     * @brief An element that performs a single optimization step.
     *        This element has its own scope and can contain child elements that run after prepare_step() but before finish_step().
     */
    class OptimizeStepElement : public LoopElement {
        public:
            OptimizeStepElement(LoopElement* owner);
            ~OptimizeStepElement() override;

            void run() override;

            /**
             * @brief Check if the current step was accepted (improved).
             */
            bool was_accepted() const {return step_accepted;}

        private:
            bool step_accepted = false;
    };
}