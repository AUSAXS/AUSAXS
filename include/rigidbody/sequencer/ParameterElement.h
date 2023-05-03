#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/LoopElementCallback.h>

#include <iostream>

namespace rigidbody {
    namespace sequencer {
        class LoopElement;
        class ParameterElement : public LoopElementCallback {
            public:
                ParameterElement(LoopElement* owner) : LoopElementCallback(owner), strategy(settings::rigidbody::parameter_generation_strategy) {}
                ParameterElement(LoopElement* owner, settings::rigidbody::ParameterGenerationStrategyChoice strategy) : LoopElementCallback(owner), strategy(strategy) {}
                virtual ~ParameterElement() = default;

                void apply() {
                    std::cout << "ParameterElement::apply()" << std::endl;
                }

            private:
                settings::rigidbody::ParameterGenerationStrategyChoice strategy;
        };
    }
}