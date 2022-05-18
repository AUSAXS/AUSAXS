#pragma once

#include <minimizer/Minimizer.h>

namespace mini {
    /**
     * @brief A wrapper class for the minimization algorithms provided by ROOT. 
     */
    class ROOT : public Minimizer {
        public:
            /**
             * @brief Default constructor.
             */
            ROOT() {}

            /**
             * @brief Constructor.
             * 
             * Construct a new minimizer using one of the algorithms provided by ROOT. 
             * 
             * @param package The package containing the algorithm.
             * @param algorithm The name of the algorithm. 
             */
            ROOT(std::string package, std::string algorithm);

            /**
             * @brief Constructor.
             * 
             * Construct a new minimizer using one of the algorithms provided by ROOT. 
             * 
             * @param package The package containing the algorithm.
             * @param algorithm The name of the algorithm. 
             */
            ROOT(std::string package, std::string algorithm, double(&f)(double*), Parameter param = Parameter());

            /**
             * @brief Perform the minimization.
             */
            Result minimize() const override;

            /**
             * @brief Add a parameter.
             */
            void add_parameter(const Parameter& param) override;

            /**
             * @brief Generate a landscape of the function.
             */
            Dataset landscape(unsigned int evals = 100) const override;

        private: 
            
    };
}