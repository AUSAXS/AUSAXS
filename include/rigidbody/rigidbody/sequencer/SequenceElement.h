#pragma once

namespace rigidbody {
    namespace sequencer {
        class SequenceElement {
            public:
                virtual void execute() = 0;
        };
    }
}