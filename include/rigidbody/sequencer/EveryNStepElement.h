#pragma once

#include <rigidbody/sequencer/LoopElement.h>
#include <utility/observer_ptr.h>

namespace rigidbody::sequencer {
    class EveryNStepElement : public LoopElement {
        public:
            EveryNStepElement(observer_ptr<LoopElement> owner, unsigned int n);
            ~EveryNStepElement() override;

            void run() override;

        private:
            int n;
            int loop_counter;
    };
}