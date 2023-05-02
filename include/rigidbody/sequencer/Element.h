#pragma once

#include <rigidbody/sequencer/SequenceElement.h>

namespace rigidbody {
    namespace sequencer {
        class Element : public SequenceElement {
            public:
                Element() = default;
                virtual ~Element() = default;

                virtual void execute() = 0;
        };
    }
}