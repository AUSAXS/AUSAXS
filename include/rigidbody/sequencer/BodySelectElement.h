#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/LoopElementCallback.h>

#include <iostream>

namespace rigidbody {
    namespace sequencer {
        class LoopElement;
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