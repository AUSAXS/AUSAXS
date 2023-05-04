#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/ParameterElementCallback.h>

namespace rigidbody {
    namespace sequencer {
        class LoopElement;
        class DecayElement : public ParameterElementCallback {
            public:
                DecayElement(LoopElement* owner);
                DecayElement(LoopElement* owner, settings::rigidbody::DecayStrategyChoice strategy);
                ~DecayElement();

                void apply();

                settings::rigidbody::DecayStrategyChoice strategy;
        };
    }
}