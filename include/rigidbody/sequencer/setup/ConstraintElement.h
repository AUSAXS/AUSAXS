#pragma once

#include <rigidbody/sequencer/GenericElement.h>

namespace rigidbody::sequencer {
    class ConstraintElement : public GenericElement {
        public:
            ConstraintElement() = default;
            ~ConstraintElement() override = default;

            void run() override;
    };
}