#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <settings/RigidBodySettings.h>

#include <memory>
#include <vector>

namespace rigidbody {
    namespace sequencer {
        /**
         * @brief A loop element is a sequence element that repeats whatever is inside it a number of times.
         */
        class LoopElement {
            public:
                LoopElement();
                LoopElement(LoopElement* owner);
                LoopElement(LoopElement* owner, unsigned int repeats);
                virtual ~LoopElement();

                virtual void execute();

                /**
                 * @brief Create a nested loop.
                 */
                LoopElement& loop(unsigned int repeats);

                /**
                 * @brief Set the parameter strategy.
                 */
                ParameterElement& parameter_strategy(settings::rigidbody::ParameterGenerationStrategyChoice strategy);

                /**
                 * @brief Set the body selection strategy.
                 */
                BodySelectElement& body_select_strategy(settings::rigidbody::BodySelectStrategyChoice strategy);

                /**
                 * @brief Set the transformation strategy.
                 */
                TransformElement& transform_strategy(settings::rigidbody::TransformationStrategyChoice strategy);

                /**
                 * @brief End the current loop.
                 */
                LoopElement& end();

                /**
                 * @brief Run an iteration of this loop. 
                 */
                void run();

                LoopElement* owner;
            protected: 
                unsigned int iterations = 1;
                std::vector<std::unique_ptr<LoopElement>> inner_loops;
                std::unique_ptr<ParameterElement> parameter_element;
                std::unique_ptr<BodySelectElement> body_select_element;
                std::unique_ptr<TransformElement> transform_element;
        };
    }
}