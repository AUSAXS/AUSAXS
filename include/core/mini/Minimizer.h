// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <mini/MiniFwd.h>
#include <mini/detail/Landscape.h>
#include <mini/detail/Result.h>

#include <functional>
#include <vector>

namespace ausaxs::mini {
    /**
     * @brief A common interface for global minimizers. 
     */
    class Minimizer {
        public:
            /**
             * @brief The function to be minimized, given as its vector of residuals r(p). The minimized quantity is chi2 = sum_i r_i(p)^2.
             */
            using residual_function = std::function<std::vector<double>(const std::vector<double>&)>;

            Minimizer();
            Minimizer(residual_function function);
            virtual ~Minimizer();

            /**
             * @brief Set the function to be minimized.
             */
            void set_function(residual_function function);

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
            [[nodiscard]] virtual mini::Landscape landscape(int bins = 100);

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
            virtual void set_max_evals(int evals);

            double tol = 1e-4;
        protected:
            std::vector<Parameter> parameters;
            mini::Landscape evaluations;
            int fevals = 0;
            int max_evals = 100;

            /**
             * @brief Evaluate the residuals at the given point.
             *        The evaluation is recorded with its chi2 unless recording has been disabled.
             */
            [[nodiscard]] std::vector<double> residuals(const std::vector<double>& params);

            /**
             * @brief Evaluate chi2, the sum of the squared residuals, at the given point.
             *        The evaluation is recorded unless recording has been disabled.
             */
            double function(const std::vector<double>& params);

            /**
             * @brief Get the sum of the squares of a residual vector.
             */
            [[nodiscard]] static double chi2(const std::vector<double>& r);

            /**
             * @brief Get the residual function of this minimizer, such that it can be handed to another minimizer.
             *        Evaluations made through it are recorded by both.
             */
            [[nodiscard]] residual_function get_recording_function();

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
            residual_function objective;
            bool record = true;

            /**
             * @brief The minimization function to be defined by subclasses. 
             *        Should be kept private, such that it can only be accessed through the common minimize() function defined here.
             */
            [[nodiscard]] virtual Result minimize_override() = 0;
    };
}