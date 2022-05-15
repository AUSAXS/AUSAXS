#pragma once

#include <utility/Axis.h>

#include <string>
#include <vector>
#include <map>
#include <functional>

/**
 * @brief A common interface for global minimizers. 
 */
class Minimizer {
    public:
        class Result {
            public:
                /**
                 * @brief Get the optimal parameters. 
                 */
                std::vector<double> parameters() const;

                /**
                 * @brief Get the errors of the parameters. 
                 */
                std::vector<double> errors() const;

            private: 
                std::map<std::string, double> params;
                std::map<std::string, double> errs;
        };

        /**
         * @brief Set the function to be minimized.
         */
        void set_function(double(&f)(double*));

        /**
         * @brief Perform the minimization.
         */
        virtual Result minimize() const = 0;

        /**
         * @brief Add a parameter.
         * 
         * @param par The name of the parameter.
         * @param guess The start value of the parameter. 
         */
        virtual void add_parameter(std::string par, double guess) = 0;

        /**
         * @brief Add a parameter with limits.
         * 
         * @param par The name of the parameter.
         * @param guess The start value of the parameter. 
         * @param limits The limits of the parameter.
         */
        virtual void add_parameter(std::string par, double guess, Limit limits) = 0;

    private:
        std::function<double(double*)> function;
        std::map<std::string, double> params;
};