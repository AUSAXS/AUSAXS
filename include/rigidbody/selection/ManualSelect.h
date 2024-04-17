#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

namespace rigidbody {
    namespace selection {
        /**
         * @brief The next body is manually selected, with the next constraint being a random one from the constraints connecting to that body.
         */
        class ManualSelect : public BodySelectStrategy {
            public: 
                ManualSelect(const RigidBody* rigidbody);
                ~ManualSelect() override;

                std::pair<unsigned int, int> next() override; ///< @copydoc BodySelectStrategy::next()

            private:
                unsigned int ibody = 0; 		// The index of the body to be transformed. 
                unsigned int iconstraint = 0; 	// The index of the constraint to be transformed.
        };
    }
}