#pragma once

#include <utility/Axis.h>
#include <utility/Dataset.h>
#include <minimizer/Utility.h>

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
            Minimizer(double(&function)(const double*));

            /**
             * @brief Constructor.
             * 
             * Initialize this minimizer with a function and dimensionality.
             */
            Minimizer(std::function<double(const double*)> function);

            /**
             * @brief Set the function to be minimized.
             */
            virtual void set_function(double(&function)(const double*));

            /**
             * @brief Set the function to be minimized.
             */
            virtual void set_function(std::function<double(const double*)> function);

            /**
             * @brief Perform the minimization.
             */
            virtual Result minimize() = 0;

            /**
             * @brief Add a parameter.
             */
            virtual void add_parameter(const Parameter& param);

            /**
             * @brief Generate a landscape of the function values. 
             *        Only valid for 1D or 2D problems.
             */
            virtual Dataset landscape(unsigned int bins) const = 0;

            /**
             * @brief Change whether the evaluations are recorded or not.
             */
            void record_evaluations(bool setting);

            double tol = 1e-6;
        protected:
            std::vector<Parameter> parameters;
            std::function<double(const double*)> function;
            std::vector<Evaluation> evaluations;

            /**
             * @brief Check if the function is set.
             */
            bool is_function_set() const noexcept;

            /**
             * @brief Check if at least one parameter has been provided.
             */
            bool is_parameter_set() const noexcept;

        private:
            std::function<double(const double*)> wrapper;
            std::function<double(const double*)> raw;
    };
}