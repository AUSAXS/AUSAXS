#pragma once

#include <rigidbody/sequencer/SequenceElement.h>
#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/sequencer/TransformElement.h>
#include <settings/RigidBodySettings.h>

#include <memory>
#include <vector>

namespace rigidbody {
    namespace sequencer {
        class LoopElement {
            public:
                LoopElement(LoopElement* owner);
                LoopElement(unsigned int iterations);
                virtual ~LoopElement() = default;

                virtual void execute();

                LoopElement& loop(unsigned int repeats);

                ParameterElement& parameter_strategy(settings::rigidbody::ParameterGenerationStrategyChoice strategy);

                BodySelectElement& body_select_strategy(settings::rigidbody::BodySelectStrategyChoice strategy);

                TransformElement& transform_strategy(settings::rigidbody::TransformationStrategyChoice strategy);

                void run();

            protected: 
                LoopElement* owner;
                unsigned int iterations = 1;
                std::vector<std::unique_ptr<LoopElement>> inner_loops;
                std::unique_ptr<ParameterElement> parameter_element;
                std::unique_ptr<BodySelectElement> body_select_element;
                std::unique_ptr<TransformElement> transform_element;
        };
    }
}