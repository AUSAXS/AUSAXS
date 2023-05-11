#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/ParameterElementCallback.h>

namespace rigidbody {
    namespace sequencer {
        class ParameterElement;
        class DecayElement : public ParameterElementCallback {
            public:
                DecayElement(rigidbody::sequencer::ParameterElement* owner);
                DecayElement(rigidbody::sequencer::ParameterElement* owner, const settings::rigidbody::DecayStrategyChoice& strategy);
                ~DecayElement() override;

                void apply();

                settings::rigidbody::DecayStrategyChoice strategy;
        };
    }
}