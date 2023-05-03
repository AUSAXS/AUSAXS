#pragma once

#include <settings/RigidBodySettings.h>

namespace rigidbody {
    namespace sequencer {
        class LoopElement;
        class ParameterElement;
        class BodySelectElement;
        class TransformElement;
        class LoopElementCallback {
            public: 
                LoopElementCallback(LoopElement* caller);

                LoopElement& loop(unsigned int repeats);

                ParameterElement& parameter_strategy(settings::rigidbody::ParameterGenerationStrategyChoice strategy);

                BodySelectElement& body_select_strategy(settings::rigidbody::BodySelectStrategyChoice strategy);

                TransformElement& transform_strategy(settings::rigidbody::TransformationStrategyChoice strategy);

                void execute();

            private: 
                LoopElement* caller;
        };
    }
}