#pragma once

#include <utility/Concepts.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>

namespace rigidbody {
    /**
     * @brief Thread-safe parameter generation strategy. The current step size scales linearly with the iteration number. 
     */
    class RotationsOnly : public ParameterGenerationStrategy {
        public: 
            /**
             * @brief Constructor.
             * 
             * @param iterations The expected number of iterations. 
             * @param length_start The start length of the generated translation vectors. 
             * @param rad_start The start angle in radians of the generated rotations. 
             */
            RotationsOnly(int iterations, double length_start, double rad_start);

            /**
             * @brief Destructor.
             */
            ~RotationsOnly() override;

            std::tuple<double, double, double> get_rotation() override;

            Vector3<double> get_translation() override;

        private: 
            double scale() const;
    };
}