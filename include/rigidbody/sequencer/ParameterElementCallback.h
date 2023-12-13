#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>

namespace rigidbody {
    namespace sequencer {
        class ParameterElementCallback : public LoopElementCallback {
            public:
                ParameterElementCallback(ParameterElement* caller);
                virtual ~ParameterElementCallback() override;

                ParameterElement& amplitude(double amplitude);

            private:
                ParameterElement* caller;
        };
    }
}