#pragma once

#include <fitter/refactored/Fitter.h>

#include <vector>

namespace fitter {
    /**
     * @brief A simple linear least-squares fitter for fitting the linear relationship y = ax+b.
     */
    class SimpleLeastSquares : public Fitter {
        public:
            /**
             * @brief Prepare a linear least-squares fit with unity errors. 
             *
             * @param data The measured data.
             * @param model The model data to be fitted.
             */
            SimpleLeastSquares(const std::vector<double> data, const std::vector<double> model);

            /**
             * @brief Prepare a linear least-squares fit.
             *
             * @param data The measured data.
             * @param model The model data to be fitted.
             * @param errors The errors on the measured data.
             */
            SimpleLeastSquares(const std::vector<double> data, const std::vector<double> model, const std::vector<double> errors);
            virtual ~SimpleLeastSquares() override = default;

            [[nodiscard]] virtual std::unique_ptr<FitResult> fit() override;
            [[nodiscard]] virtual double fit_chi2_only() override;
            [[nodiscard]] unsigned int dof() const override;
            [[nodiscard]] unsigned int size() const override;

        private:
            std::vector<double> data, model, inv_sigma2;

            /**
             * @brief Perform a linear least-squares fit and calculate @a only the fitted parameters.
             *
             * @return The fitted parameters (a, b, a_err^2, b_err^2) for the equation y = ax+b.
             */
            [[nodiscard]] std::array<double, 4> fit_params_only();

            /**
             * @brief Calculate chi2.
             */
            [[nodiscard]] double chi2(const std::vector<double>& params) override;
    };
}