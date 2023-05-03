#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/LoopElementCallback.h>

#include <iostream>

namespace rigidbody {
    namespace sequencer {
        class LoopElement;
        class TransformElement : public LoopElementCallback {
            public:
                TransformElement(LoopElement* owner) : LoopElementCallback(owner), strategy(settings::rigidbody::transform_strategy) {}
                TransformElement(LoopElement* owner, settings::rigidbody::TransformationStrategyChoice strategy) : LoopElementCallback(owner), strategy(strategy) {}
                ~TransformElement() = default;

                void apply() {
                    std::cout << "TransformElement::apply()" << std::endl;
                }

            private:
                settings::rigidbody::TransformationStrategyChoice strategy;
        };
    }
}