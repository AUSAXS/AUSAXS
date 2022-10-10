#pragma once

#include <utility/Axis.h>
#include <utility/Dataset2D.h>
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
             * @brief Destructor.
             */
            virtual ~Minimizer() = default;

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
            Result minimize();

            /**
             * @brief Add a parameter.
             */
            virtual void add_parameter(const Parameter& param);

            /**
             * @brief Remove any set parameters.
             */
            void clear_parameters() noexcept;

            /**
             * @brief Generate a landscape of the function values. 
             *        Only valid for 1D or 2D problems.
             */
            virtual Dataset2D landscape(unsigned int bins = 100) = 0;

            /**
             * @brief Get the evaluated points. 
             */
            virtual Dataset2D get_evaluated_points() const = 0;

            /**
             * @brief Check if this minimizer has been initialized.
             */
            bool empty() const noexcept;

            /**
             * @brief Change whether the evaluations are recorded or not.
             */
            void record_evaluations(bool setting);

            double tol = 1e-4;
        protected:
            std::vector<Parameter> parameters;
            std::function<double(const double*)> function;
            std::vector<Evaluation> evaluations;
            unsigned int fevals = 0;

            /**
             * @brief Clear the evaluated points.
             */
            void clear_evaluated_points() noexcept;

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

            /**
             * @brief The minimization function to be defined by subclasses. 
             *        Should be kept private, such that it can only be accessed through the common minimize() function defined here.
             */
            virtual Result minimize_override() = 0;
    };
}