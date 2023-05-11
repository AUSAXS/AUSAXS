#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/LoopElementCallback.h>

namespace rigidbody {
    namespace sequencer {
        class LoopElement;
        class DecayElement;
        class ParameterElement : public LoopElementCallback {
            public:
                ParameterElement(LoopElement* owner);
                ParameterElement(LoopElement* owner, settings::rigidbody::ParameterGenerationStrategyChoice strategy);
                ~ParameterElement() override;

                void apply();

                ParameterElement& amplitude(double amplitude);

                DecayElement& decay_strategy(const settings::rigidbody::DecayStrategyChoice& strategy);

                LoopElement* owner;
                settings::rigidbody::ParameterGenerationStrategyChoice strategy;
            private: 
                std::unique_ptr<DecayElement> decay_element;

                void initialize();
        };
    }
}