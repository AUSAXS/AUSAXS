#pragma once

#include <fitter/Fitter.h>

#include <vector>

namespace fitter {
    /**
     * @brief A simple linear least-squares fitter for fitting the linear relationship y = ax+b.
     */
    class LinearFitter : public Fitter {
        public:
            virtual ~LinearFitter() override = default;

            /**
             * @brief Prepare a linear least-squares fit with unity errors. 
             *
             * @param data The measured data.
             * @param model The model data to be fitted.
             */
            LinearFitter(const std::vector<double> data, const std::vector<double> model);

            /**
             * @brief Prepare a linear least-squares fit.
             *
             * @param data The measured data.
             * @param model The model data to be fitted.
             * @param errors The errors on the measured data.
             */
            LinearFitter(const std::vector<double> data, const std::vector<double> model, const std::vector<double> errors);

            [[nodiscard]] std::unique_ptr<FitResult> fit() override;
            [[nodiscard]] double fit_chi2_only() override;
            [[nodiscard]] unsigned int dof() const override;
            [[nodiscard]] unsigned int size() const override;
            [[nodiscard]] std::vector<double> get_residuals(const std::vector<double>& params) override;

        private:
            std::vector<double> data, model, inv_sigma;

            /**
             * @brief Perform a linear least-squares fit and calculate @a only the fitted parameters.
             *
             * @return The fitted parameters (a, b, a_err^2, b_err^2) for the equation y = ax+b.
             */
            [[nodiscard]] std::vector<double> fit_params_only();
    };
}