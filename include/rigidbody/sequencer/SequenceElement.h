#pragma once

namespace rigidbody {
    namespace sequencer {
        class SequenceElement {
            public:
                SequenceElement() = default;
                virtual ~SequenceElement() = default;

                virtual void execute() = 0;
        };
    }
}