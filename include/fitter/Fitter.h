#pragma once

#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <utility>
#include <memory>

#include <dataset/Multiset.h>
#include <dataset/Dataset.h>
#include <fitter/Fit.h>

namespace fitter {
    class Fitter {
        public:
            virtual ~Fitter() {}

            /**
             * @brief Perform a fit and return a fit object containing various information. 
             */
            [[nodiscard]] virtual std::shared_ptr<Fit> fit() = 0;

            /**
             * @brief Perform a fit and return the minimum function value.
             */
            [[nodiscard]] virtual double fit_only() = 0;

            /**
             * @brief Make a plot of the fit. 
             * 
             * @return A vector of TGraphs {Interpolated points, Optimal line, Measured points with uncertainties}
             */
            [[nodiscard]] virtual Fit::Plots plot() = 0;

            /**
             * @brief Make a residual plot of the fit.
             * 
             * @return A TGraphErrors with the residuals and their uncertainties. 
             */
            [[nodiscard]] virtual SimpleDataset plot_residuals() = 0;

            [[nodiscard]] virtual std::shared_ptr<Fit> get_fit() const = 0;

            [[nodiscard]] virtual unsigned int dof() const = 0;

            [[nodiscard]] virtual double chi2(const std::vector<double>& params) = 0;
    };
}