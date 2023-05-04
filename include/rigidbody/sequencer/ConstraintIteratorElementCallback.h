#pragma once

#include <rigidbody/sequencer/LoopElementCallback.h>

namespace rigidbody {
    namespace sequencer {
        class ConstraintIteratorElement;
        class ConstraintIteratorElementCallback : public LoopElementCallback {
            public:
                ConstraintIteratorElementCallback(ConstraintIteratorElement* caller);

                BodySelectElement& body_select_strategy(settings::rigidbody::BodySelectStrategyChoice strategy) = delete;

            private:
                ConstraintIteratorElement* caller;
        };
    }
}