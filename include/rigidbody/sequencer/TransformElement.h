#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <settings/RigidBodySettings.h>

namespace rigidbody {
    namespace sequencer {
        class TransformElement : public LoopElementCallback {
            public:
                TransformElement(LoopElement* owner);
                TransformElement(LoopElement* owner, settings::rigidbody::TransformationStrategyChoice strategy);
                ~TransformElement();

                void apply();

                settings::rigidbody::TransformationStrategyChoice strategy;
        };
    }
}