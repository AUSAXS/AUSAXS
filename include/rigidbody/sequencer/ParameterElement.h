#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/DecayElement.h>

namespace rigidbody {
    namespace sequencer {
        class LoopElement;
        class ParameterElement : public LoopElementCallback {
            public:
                ParameterElement(LoopElement* owner);
                ParameterElement(LoopElement* owner, settings::rigidbody::ParameterGenerationStrategyChoice strategy);
                ~ParameterElement();

                void apply();

                ParameterElement& amplitude(double amplitude);

                DecayElement& decay_strategy(settings::rigidbody::DecayStrategyChoice strategy);

                LoopElement* owner;
                settings::rigidbody::ParameterGenerationStrategyChoice strategy;
            private: 
                std::unique_ptr<DecayElement> decay_element = std::make_unique<DecayElement>();                
        };
    }
}