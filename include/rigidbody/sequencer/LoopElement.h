#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/GenericElement.h>
#include <fitter/FitterFwd.h>

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

                virtual std::shared_ptr<fitter::Fit> execute();

                /**
                 * @brief Create a nested loop.
                 */
                LoopElement& loop(unsigned int repeats);

                /**
                 * @brief Set the parameter strategy.
                 */
                ParameterElement& parameter_strategy(std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy);

                /**
                 * @brief Set the body selection strategy.
                 */
                BodySelectElement& body_select_strategy(std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy);

                /**
                 * @brief Set the transformation strategy.
                 */
                TransformElement& transform_strategy(std::unique_ptr<rigidbody::transform::TransformStrategy> strategy);

                /**
                 * @brief Perform a single optimization step.
                 */
                LoopElement& optimize();

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
                std::vector<std::unique_ptr<GenericElement>> elements;
        };
    }
}