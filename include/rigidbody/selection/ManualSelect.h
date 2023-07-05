#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

namespace rigidbody {
    namespace selection {
        /**
         * @brief Thread-safe body selection strategy. The next body is manually selected.
         */
        class ManualSelect : public BodySelectStrategy {
            public: 
                /**
                 * @brief Constructor.
                 */
                ManualSelect(const RigidBody* rigidbody);

                /**
                 * @brief Destructor.
                 */
                ~ManualSelect() override;

                /**
                 * @brief Get the index of the next body to be transformed. 
                 */
                std::pair<unsigned int, unsigned int> next() override;

            private:
                unsigned int ibody = 0; 		// The index of the body to be transformed. 
                unsigned int iconstraint = 0; 	// The index of the constraint to be transformed.
        };
    }
}