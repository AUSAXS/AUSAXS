#pragma once

#include <minimizer/Minimizer.h>

class ScanMinimizer : public Minimizer {
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