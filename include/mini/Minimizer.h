#pragma once

#include <utility/Axis.h>
#include <dataset/Dataset2D.h>
#include <mini/Utility.h>

#include <string>
#include <vector>
#include <map>
#include <functional>

namespace mini {
    enum class type {
        BFGS,
        DLIB_GLOBAL,
        GOLDEN,
        MINIMUM_EXPLORER,
        SCAN,
        LIMITED_SCAN
    };

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
            Minimizer(double(&function)(std::vector<double>));

            /**
             * @brief Constructor.
             * 
             * Initialize this minimizer with a function and dimensionality.
             */
            Minimizer(std::function<double(std::vector<double>)> function);

            /**
             * @brief Destructor.
             */
            virtual ~Minimizer() = default;

            /**
             * @brief Set the function to be minimized.
             */
            virtual void set_function(double(&function)(std::vector<double>));

            /**
             * @brief Set the function to be minimized.
             */
            virtual void set_function(std::function<double(std::vector<double>)> function);

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
            virtual mini::Landscape landscape(unsigned int bins = 100);

            /**
             * @brief Get the evaluated points. 
             */
            mini::Landscape get_evaluated_points() const;

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
            std::function<double(std::vector<double>)> function = [](std::vector<double>){throw except::unexpected("Minimizer::function: Function was not initialized."); return 0.0;};
            mini::Landscape evaluations;
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
            std::function<double(std::vector<double>)> wrapper;
            std::function<double(std::vector<double>)> raw;

            /**
             * @brief The minimization function to be defined by subclasses. 
             *        Should be kept private, such that it can only be accessed through the common minimize() function defined here.
             */
            virtual Result minimize_override() = 0;
    };
}