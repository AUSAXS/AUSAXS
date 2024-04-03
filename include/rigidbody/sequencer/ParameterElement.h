#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/GenericElement.h>

namespace rigidbody {
    namespace sequencer {
        class ParameterElement : public LoopElementCallback, public GenericElement {
            public:
                ParameterElement(LoopElement* owner);
                ParameterElement(LoopElement* owner, settings::rigidbody::ParameterGenerationStrategyChoice strategy);
                ~ParameterElement() override;

                void run() override;

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