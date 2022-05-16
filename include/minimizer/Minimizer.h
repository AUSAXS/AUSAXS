#pragma once

#include <utility/Axis.h>
#include <utility/Dataset.h>

#include <string>
#include <vector>
#include <map>
#include <functional>

namespace mini {
    /**
     * @brief A representation of a parameter.
     */
    struct Parameter {
        Parameter(std::string name, double guess, const Limit& bounds) : name(name), guess(guess), bounds(bounds) {}

        bool is_bounded() const {return bounds.empty();}

        std::string name;
        double guess;
        Limit bounds;
    };

    /**
     * @brief A common interface for global minimizers. 
     */
    class Minimizer {
        public:
            struct Result {
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

                std::map<std::string, double> parameters;
                std::map<std::string, double> errors;
                double fval;
            };

            struct Evaluation {
                Evaluation(std::vector<double> vals, double fval) : vals(vals), fval(fval) {}

                std::vector<double> vals;
                double fval;
            };

            /**
             * @brief Default constructor.
             */
            Minimizer() {}

            /**
             * @brief Constructor.
             * 
             * Initialize this minimizer with a function and dimensionality.
             */
            Minimizer(double(&function)(double*), unsigned int dimensionality);

            /**
             * @brief Constructor.
             * 
             * Initialize this minimizer with a function and dimensionality.
             */
            Minimizer(std::function<double(double*)> function, unsigned int dimensionality);

            /**
             * @brief Set the function to be minimized.
             */
            virtual void set_function(double(&f)(double*), unsigned int dimensionality = 0);

            /**
             * @brief Set the function to be minimized.
             */
            virtual void set_function(std::function<double(double*)> function, unsigned int dimensionality = 0);

            /**
             * @brief Perform the minimization.
             */
            virtual Result minimize() const = 0;

            /**
             * @brief Add a parameter.
             * 
             * @param param The parameter. 
             */
            virtual void add_parameter(const Parameter& param) = 0;

            /**
             * @brief Change whether the evaluations are recorded or not.
             */
            void record_evaluations(bool setting);

            double tol = 1e-6;
        protected:
            std::function<double(double*)> function;
            unsigned int dimensionality;
            std::vector<Evaluation> evaluations;

        private:
            std::function<double(double*)> wrapper;
            std::function<double(double*)> raw;

    };
}