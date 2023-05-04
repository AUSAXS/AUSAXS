#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/LoopElementCallback.h>

#include <iostream>

namespace rigidbody {
    namespace sequencer {
        class LoopElement;
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