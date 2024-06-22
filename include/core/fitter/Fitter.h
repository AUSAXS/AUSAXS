#pragma once

#include <fitter/FitterFwd.h>
#include <fitter/FitResult.h>
#include <dataset/SimpleDataset.h>
#include <mini/Minimizer.h>

#include <vector>
#include <memory>

namespace fitter {
    class Fitter {
        public:
            virtual ~Fitter() = default;

            /**
             * @brief Perform a fit and return a fit object containing various information. 
             */
            [[nodiscard]] virtual std::shared_ptr<FitResult> fit() = 0;

            /**
             * @brief Perform a fit and return the minimum function value.
             */
            [[nodiscard]] virtual double fit_chi2_only() = 0;

            /**
             * @brief Plot the fitted function and the measured data points.
             */
            [[nodiscard]] virtual FitResult::FitInfo plot() = 0;

            /**
             * @brief Make a residual plot of the fit.
             */
            [[nodiscard]] virtual SimpleDataset plot_residuals() = 0;

            [[nodiscard]] virtual std::shared_ptr<FitResult> get_fit() const = 0;

            /**
             * @brief Get the number of degrees of freedom.
             */
            [[nodiscard]] virtual unsigned int dof() const = 0;

            /**
             * @brief Get the number of degrees of freedom.
             */
            [[nodiscard]] unsigned int degrees_of_freedom() const {return dof();}

            /**
             * @brief Get the total number of data points. 
             */
            [[nodiscard]] virtual unsigned int size() const = 0;

            /**
             * @brief Evaluate the chi2 function for the given parameters.
             */
            [[nodiscard]] virtual double chi2(const std::vector<double>& params) = 0;
    };
}