#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>

namespace rigidbody::sequencer {
    class ConstraintIteratorElementCallback : public LoopElementCallback {
        public:
            ConstraintIteratorElementCallback(ConstraintIteratorElement* caller);

            BodySelectElement& body_select_strategy(settings::rigidbody::BodySelectStrategyChoice strategy) = delete;

        private:
            ConstraintIteratorElement* caller;
    };
}