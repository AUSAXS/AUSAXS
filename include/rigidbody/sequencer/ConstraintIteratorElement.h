#pragma once

#include <rigidbody/sequencer/LoopElement.h>

namespace rigidbody::sequencer {
    /**
     * @brief Iterates over all constraints in the body.
     */
    class ConstraintIteratorElement : public LoopElement {
        public:
            using LoopElement::LoopElement;
            ~ConstraintIteratorElement() = default;

            void execute() override;
    };
}