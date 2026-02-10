// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/elements/LoopElementCallback.h>

namespace ausaxs::rigidbody::sequencer {
    class ConstraintIteratorElementCallback : public LoopElementCallback {
        public:
            ConstraintIteratorElementCallback(ConstraintIteratorElement* caller);

            BodySelectElement& body_select_strategy(settings::rigidbody::BodySelectStrategyChoice strategy) = delete;

        private:
            ConstraintIteratorElement* caller;
    };
}