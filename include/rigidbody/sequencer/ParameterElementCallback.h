#pragma once

#include <rigidbody/sequencer/LoopElementCallback.h>

namespace rigidbody {
    namespace sequencer {
        class ParameterElement;
        class ParameterElementCallback : public LoopElementCallback {
            public:
                ParameterElementCallback(ParameterElement* caller);

                ParameterElement& amplitude(double amplitude);

            private:
                ParameterElement* caller;
        };
    }
}