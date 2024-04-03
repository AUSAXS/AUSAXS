#pragma once

#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>

namespace rigidbody::sequencer {
    class OptimizeStepElement : public LoopElementCallback, public GenericElement {
        public:
            OptimizeStepElement(LoopElement* owner);
            ~OptimizeStepElement() override;

            void run() override;
    };
}