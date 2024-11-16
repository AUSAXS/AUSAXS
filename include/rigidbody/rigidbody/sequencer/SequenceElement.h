#pragma once

namespace ausaxs::rigidbody::sequencer {
    class SequenceElement {
        public:
            virtual void execute() = 0;
    };
}