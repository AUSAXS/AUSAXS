#pragma once

#include <utility/Axis.h>

#include <string>
#include <vector>
#include <map>
#include <functional>

namespace mini {
    /**
     * @brief A common interface for global minimizers. 
     */
    class Minimizer {
        public:
            class Result {
                public:
                    /**
                     * @brief Default constructor.
                     */
                    Result() {}

                    /**
                     * @brief Construct a Result. 
                     * 
                     * @param params The parameter values. 
                     * @param errors The parameter errors. 
                     * @param function_val The function value.
                     */
                    Result(const std::map<std::string, double>& params, const std::map<std::string, double>& errors, double function_val);

                    /**
                     * @brief Construct a Result. 
                     * 
                     * @param values The parameter values. 
                     * @param errors The parameter errors. 
                     * @param names The parameter names.
                     * @param function_val The function value.
                     */
                    Result(const std::vector<double>& values, const std::vector<double>& errors, const std::vector<std::string>& names, double function_val);

                    /**
                     * @brief Get the optimal parameters. 
                     */
                    std::map<std::string, double> parameters() const;

                    /**
                     * @brief Get the errors of the parameters. 
                     */
                    std::map<std::string, double> errors() const;

                private: 
                    std::map<std::string, double> params;
                    std::map<std::string, double> errs;
                    double fval;
            };

            /**
             * @brief Default constructor.
             */
            Minimizer() {}

            /**
             * @brief Constructor.
             * 
             * Initialize this minimizer with a function.
             */
            Minimizer(double(&f)(double*));

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

        protected:
            double tol = 1e-6;
            std::function<double(double*)> function;
            std::vector<double> params;
            std::vector<std::string> param_names;
    };
}