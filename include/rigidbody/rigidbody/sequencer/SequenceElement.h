#pragma once

namespace rigidbody::sequencer {
    class SequenceElement {
        public:
            virtual void execute() = 0;
    };
}