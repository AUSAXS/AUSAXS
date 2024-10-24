#pragma once

#include <fitter/FitterFwd.h>
#include <fitter/FitResult.h>
#include <mini/MiniFwd.h>

#include <vector>
#include <memory>

namespace ausaxs::fitter {
    class Fitter {
        public:
            virtual ~Fitter() = default;

            /**
             * @brief Perform a fit and return a fit object containing various information. 
             */
            [[nodiscard]] virtual std::unique_ptr<FitResult> fit() = 0;

            /**
             * @brief Perform a fit and return the minimum function value.
             */
            [[nodiscard]] double fit_chi2_only();

            /**
             * @brief Perform a fit and return the optimal parameters.
             */
            [[nodiscard]] virtual std::vector<double> fit_params_only() = 0;

            /**
             * @brief Get the number of degrees of freedom.
             */
            [[nodiscard]] virtual unsigned int dof() const = 0;
            [[nodiscard]] unsigned int degrees_of_freedom() const {return dof();} //< @copydoc dof

            /**
             * @brief Get the total number of data points. 
             */
            [[nodiscard]] virtual unsigned int size() const = 0;

            /**
             * @brief Evaluate the chi2 for the given parameters.
             */
            [[nodiscard]] virtual double chi2(const std::vector<double>& params);

            /**
             * @brief Get the residuals for the given parameters.
             */
            [[nodiscard]] virtual std::vector<double> get_residuals(const std::vector<double>& params) = 0;

            /**
             * @brief Set the minimization algorithm to use.
             */
            void set_algorithm(mini::algorithm t) {algorithm = t;}

        protected:
            mini::algorithm algorithm = mini::algorithm::DEFAULT;
    };
}