#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>

namespace rigidbody {
    namespace sequencer {
        class BodySelectElement : public LoopElementCallback, public GenericElement {
            public:
                BodySelectElement(LoopElement* owner);
                BodySelectElement(LoopElement* owner, settings::rigidbody::BodySelectStrategyChoice strategy);
                ~BodySelectElement();

                void run() override;

                settings::rigidbody::BodySelectStrategyChoice strategy;
        };
    }
}