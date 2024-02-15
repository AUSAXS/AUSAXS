#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <settings/RigidBodySettings.h>

namespace rigidbody {
    namespace sequencer {
        class BodySelectElement : public LoopElementCallback {
            public:
                BodySelectElement(LoopElement* owner);
                BodySelectElement(LoopElement* owner, settings::rigidbody::BodySelectStrategyChoice strategy);
                ~BodySelectElement();

                void apply();

                settings::rigidbody::BodySelectStrategyChoice strategy;
        };
    }
}