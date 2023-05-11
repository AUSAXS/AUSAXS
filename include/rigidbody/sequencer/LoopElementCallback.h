#pragma once

#include <settings/RigidBodySettings.h>

namespace rigidbody {
    namespace sequencer {
        class LoopElement;
        class ParameterElement;
        class BodySelectElement;
        class TransformElement;

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
                 * @brief Execute this entire sequence. This is the last method that should be called.
                 */
                void execute();

            private: 
                LoopElement* caller;
        };
    }
}