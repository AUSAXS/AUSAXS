#pragma once

#include <rigidbody/sequencer/LoopElementCallback.h>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Iterates over all constraints in the body.
     */
    class ConstraintIteratorElement : public LoopElementCallback {
        public:
            using LoopElementCallback::LoopElementCallback;
            ~ConstraintIteratorElement() = default;
    };
}