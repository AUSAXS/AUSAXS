#pragma once

#include <rigidbody/sequencer/SequenceElement.h>
#include <settings/RigidBodySettings.h>

namespace rigidbody {
    namespace sequencer {
        class Parameters : public SequenceElement {
            public:
                Parameters() = default;
                Parameters(settings::rigidbody::ParameterGenerationStrategyChoice strategy) : strategy(strategy) {}
                virtual ~Parameters() = default;
                void execute();

            private: 
                settings::rigidbody::ParameterGenerationStrategyChoice strategy;
        };
    }
}