#pragma once

#include <utility/Axis.h>
#include <utility/Dataset.h>

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
             * @brief Generate a landscape of the function values. 
             *        Only valid for 1D or 2D problems.
             */
            virtual Dataset landscape() const = 0;

            /**
             * @brief Change whether the evaluations are recorded or not.
             */
            void record_evaluations(bool setting);

            double tol = 1e-6;
        protected:
            std::vector<Parameter> parameters;
            std::function<double(double*)> function;
            unsigned int dimensionality;
            std::vector<Evaluation> evaluations;

        private:
            std::function<double(double*)> wrapper;
            std::function<double(double*)> raw;

    };
}