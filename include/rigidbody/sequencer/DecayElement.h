#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/ParameterElementCallback.h>
#include <settings/RigidBodySettings.h>

namespace rigidbody {
    namespace sequencer {
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