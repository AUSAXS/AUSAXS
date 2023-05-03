#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/LoopElementCallback.h>

#include <iostream>

namespace rigidbody {
    namespace sequencer {
        class LoopElement;
        class BodySelectElement : public LoopElementCallback {
            public:
                BodySelectElement(LoopElement* owner) : LoopElementCallback(owner), strategy(settings::rigidbody::body_select_strategy) {}
                BodySelectElement(LoopElement* owner, settings::rigidbody::BodySelectStrategyChoice strategy) : LoopElementCallback(owner), strategy(strategy) {}
                virtual ~BodySelectElement() = default;

                void apply() {
                    std::cout << "BodySelectElement::apply()" << std::endl;
                }

            private:
                settings::rigidbody::BodySelectStrategyChoice strategy;
        };
    }
}