#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/RigidbodyFwd.h>

namespace rigidbody {
    namespace sequencer {
        /**
         * @brief A callback class for the LoopElement class.
         * 
         * This class is used to provide a callback interface for the LoopElement class, giving other elements access
         * to the basic loop methods.
         */
        class LoopElementCallback {
            public: 
                /**
                 * @brief Create a callback interface.
                 */
                LoopElementCallback(LoopElement* caller);

                virtual ~LoopElementCallback();

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
                 * @brief End the current loop.
                 */
                LoopElement& end();

                /**
                 * @brief Execute this entire sequence. This is the last method that should be called.
                 */
                void execute();

                LoopElement* owner;
        };
    }
}