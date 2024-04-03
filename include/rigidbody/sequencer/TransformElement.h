#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>

namespace rigidbody {
    namespace sequencer {
        class TransformElement : public LoopElementCallback, public GenericElement {
            public:
                TransformElement(LoopElement* owner);
                TransformElement(LoopElement* owner, settings::rigidbody::TransformationStrategyChoice strategy);
                ~TransformElement();

                void run() override;

                settings::rigidbody::TransformationStrategyChoice strategy;
        };
    }
}