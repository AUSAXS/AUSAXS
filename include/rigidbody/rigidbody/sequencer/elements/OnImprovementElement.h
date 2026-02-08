// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/LoopElement.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    class OptimizeStepElement;

    /**
     * @brief A conditional element that runs its children only if the optimization step was accepted.
     *        This must be a child of OptimizeStepElement.
     */
    class OnImprovementElement : public LoopElement {
        public:
            OnImprovementElement(observer_ptr<OptimizeStepElement> owner);
            ~OnImprovementElement() override;

            void run() override;

        private:
            observer_ptr<OptimizeStepElement> owner;
    };
}