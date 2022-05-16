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
            ROOT(std::string package, std::string algorithm, double(&f)(double*), std::map<std::string, double> params = {});

            /**
             * @brief Perform the minimization.
             */
            Result minimize() const override;

            /**
             * @brief Add a parameter.
             * 
             * @param par The name of the parameter.
             * @param guess The start value of the parameter. 
             */
            void add_parameter(std::string par, double guess);

            /**
             * @brief Add a parameter with limits.
             * 
             * @param par The name of the parameter.
             * @param guess The start value of the parameter. 
             * @param limits The limits of the parameter.
             */
            void add_parameter(std::string par, double guess, Limit limits);
    };
}