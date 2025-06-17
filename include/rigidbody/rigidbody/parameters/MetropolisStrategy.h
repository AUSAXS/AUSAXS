#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>

namespace ausaxs::rigidbody::sampling {
    /**
     * @brief A Metropolis sampling strategy for rigid body parameters.
     * 
     * This strategy uses the Metropolis algorithm to sample rigid body parameters.
     */
    class MetropolisStrategy : public parameter::ParameterGenerationStrategy {
        public:
            MetropolisStrategy(observer_ptr<const Rigidbody> molecule, unsigned int iterations);
            ~MetropolisStrategy() override;

            void set_pdf(std::function<double(double)> pdf) {this->pdf = std::move(pdf);}
            bool accept(double x_old, double x_new);

        private:
            std::function<double(double)> pdf;
    };
}