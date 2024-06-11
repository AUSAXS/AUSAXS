#pragma once

#include <mini/MiniFwd.h>
#include <mini/detail/Landscape.h>
#include <mini/detail/Result.h>

#include <vector>
#include <functional>

namespace mini {
    enum class type {
        GOLDEN,
        MINIMUM_EXPLORER,
        SCAN,
        LIMITED_SCAN,
        #if defined(DLIB_AVAILABLE)
            DLIB_GLOBAL,
            BFGS,
            DEFAULT=BFGS
        #else
            DEFAULT=GOLDEN
        #endif
    };

    /**
     * @brief A common interface for global minimizers. 
     */
    class Minimizer {
        public:
            /**
             * @brief Default constructor.
             */
            Minimizer();

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
            virtual ~Minimizer();

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
            [[nodiscard]] Result minimize();

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
            [[nodiscard]] virtual mini::Landscape landscape(unsigned int bins = 100);

            /**
             * @brief Get the evaluated points. 
             */
            [[nodiscard]] mini::Landscape get_evaluated_points() const;

            /**
             * @brief Check if this minimizer has been initialized.
             */
            [[nodiscard]] bool empty() const noexcept;

            /**
             * @brief Change whether the evaluations are recorded or not.
             */
            void record_evaluations(bool setting);

            /**
             * @brief Set the maximum number of evaluations.
             *        Note that this is not supported by all minimizers, in which case it will be ignored.
             */
            virtual void set_max_evals(unsigned int evals);

            double tol = 1e-4;
        protected:
            std::vector<Parameter> parameters;
            std::function<double(std::vector<double>)> function = [](std::vector<double>){throw std::runtime_error("Minimizer::function: Function was not initialized."); return 0.0;};
            mini::Landscape evaluations;
            unsigned int fevals = 0;
            unsigned int max_evals = 100;

            /**
             * @brief Clear the evaluated points.
             */
            void clear_evaluated_points() noexcept;

            /**
             * @brief Check if the function is set.
             */
            [[nodiscard]] bool is_function_set() const noexcept;

            /**
             * @brief Check if at least one parameter has been provided.
             */
            [[nodiscard]] bool is_parameter_set() const noexcept;

        private:
            std::function<double(std::vector<double>)> wrapper;
            std::function<double(std::vector<double>)> raw;

            /**
             * @brief The minimization function to be defined by subclasses. 
             *        Should be kept private, such that it can only be accessed through the common minimize() function defined here.
             */
            [[nodiscard]] virtual Result minimize_override() = 0;
    };
}