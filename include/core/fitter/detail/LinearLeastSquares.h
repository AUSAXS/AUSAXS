// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <fitter/Fitter.h>

#include <vector>

namespace ausaxs::fitter::detail {
    /**
     * @brief A simple linear least-squares fitter for fitting the linear relationship y = ax+b.
     */
    class LinearLeastSquares : public Fitter {
        protected:
            LinearLeastSquares() = default;

        public:
            virtual ~LinearLeastSquares() override = default;

            /**
             * @brief Prepare a linear least-squares fit with unity errors. 
             *
             * @param x The independent variable. When fitting a model to a measurement, this is the model curve.
             * @param y The dependent variable. When fitting a model to a measurement, this is the measured data.
             */
            LinearLeastSquares(const std::vector<double>& x, const std::vector<double>& y);

            /**
             * @brief Prepare a linear least-squares fit.
             *
             * @param x The independent variable. When fitting a model to a measurement, this is the model curve.
             * @param y The dependent variable. When fitting a model to a measurement, this is the measured data.
             * @param yerr The errors on @a y.
             */
            LinearLeastSquares(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& yerr);

            [[nodiscard]] std::unique_ptr<FitResult> fit() override;
            [[nodiscard]] unsigned int dof() const override;
            [[nodiscard]] unsigned int size() const override;
            [[nodiscard]] std::vector<double> get_residuals(const std::vector<double>& params) override;
            [[nodiscard]] std::vector<double> get_residuals() {return get_residuals(fit_params_only());}

            /**
             * @brief Get the model curve for the given parameters.
             */
            [[nodiscard]] std::vector<double> get_model_curve(const std::vector<double>& params);

            /**
             * @brief Fit and get the model curve.
             */
            [[nodiscard]] std::vector<double> get_model_curve();

        protected:
            std::vector<double> x, y, inv_sigma;

            /**
             * @brief Perform a linear least-squares fit and calculate @a only the fitted parameters.
             *
             * @return The fitted parameters (a, b, a_err^2, b_err^2) for the equation y = ax+b.
             */
            [[nodiscard]] std::vector<double> fit_params_only() override;
    };
}