#pragma once

namespace ausaxs::rigidbody::sequencer {
    class GenericElement {
        public:
            GenericElement() = default;
            virtual ~GenericElement() = default;

            virtual void run() = 0;
    };
}