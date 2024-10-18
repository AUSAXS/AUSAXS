#pragma once

#include <fitter/Fitter.h>
#include <dataset/DatasetFwd.h>

#include <vector>

namespace fitter {
    /**
     * @brief A simple linear least-squares fitter for fitting the linear relationship y = ax+b.
     */
    class LinearFitter : public Fitter {
        public:
            virtual ~LinearFitter() override = default;

            /**
             * @brief Prepare a fit of the measured values in @a input to the model described by @a h.
             */
            LinearFitter(const SimpleDataset& data, std::unique_ptr<hist::DistanceHistogram> model);

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

			/**
			 * @brief Set the scattering histogram used for the fit. 
			 */
			void set_model(std::unique_ptr<hist::DistanceHistogram> model);

        private:
            std::vector<double> data, model, inv_sigma;

            /**
             * @brief Perform a linear least-squares fit and calculate @a only the fitted parameters.
             *
             * @return The fitted parameters (a, b, a_err^2, b_err^2) for the equation y = ax+b.
             */
            [[nodiscard]] std::vector<double> fit_params_only() override;
    };
}