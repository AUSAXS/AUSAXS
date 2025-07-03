// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>

namespace ausaxs::rigidbody::sequencer {
    class OptimizeStepElement : public LoopElementCallback, public GenericElement {
        public:
            OptimizeStepElement(LoopElement* owner);
            ~OptimizeStepElement() override;

            void run() override;

        private:
            std::vector<std::unique_ptr<GenericElement>> elements;
    };
}