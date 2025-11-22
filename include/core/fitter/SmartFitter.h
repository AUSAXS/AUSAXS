// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <fitter/LinearFitter.h>
#include <dataset/SimpleDataset.h>
#include <mini/detail/Parameter.h>
#include <hist/HistFwd.h>

namespace ausaxs::fitter {
    /**
     * @brief A smart fitter automatically fitting the parameters enabled in the settings. 
     *
     * Supported parameters: linear (always), solvent shell density, excluded volume, solvent density
     */
    class SmartFitter : public Fitter {
        public:
            struct EnabledFitParameters {
                bool hydration, excluded_volume, solvent_density, atomic_debye_waller, exv_debye_waller;
                unsigned int get_enabled_pars_count() const;
                void apply_pars(const std::vector<double>& params, observer_ptr<hist::DistanceHistogram> model);
                void validate_model(observer_ptr<hist::DistanceHistogram> h);
                static EnabledFitParameters initialize_from_settings();
            } enabled_fit_parameters;

            virtual ~SmartFitter() override;
            SmartFitter(SmartFitter&&) noexcept;
            SmartFitter& operator=(SmartFitter&&) noexcept;

            /**
             * @brief Prepare a fit of the measured values in @a input to a model to be defined later. 
             */
            SmartFitter(const SimpleDataset& data);

            /**
             * @brief Prepare a fit of the measured values in @a input to the model described by @a h.
             */
            SmartFitter(const SimpleDataset& data, std::unique_ptr<hist::DistanceHistogram> h);

            [[nodiscard]] virtual std::unique_ptr<FitResult> fit() override;
            [[nodiscard]] unsigned int dof() const override;
            [[nodiscard]] unsigned int size() const override;

            /**
             * @brief Set the guess values for the fit. 
             */
            void set_guess(std::vector<mini::Parameter>&& guess);

			/**
			 * @brief Set the scattering histogram used for the fit. 
			 */
			void set_model(std::unique_ptr<hist::DistanceHistogram> h);

			/**
			 * @brief Get the scattering histogram used for the fit. 
			 */
			[[nodiscard]] observer_ptr<hist::DistanceHistogram> get_model();

            /**
             * @brief Get the model curve for the given parameters.
             */
            [[nodiscard]] std::vector<double> get_model_curve(const std::vector<double>& params);

            /**
             * @brief Get the dataset being fitted. 
             */
            SimpleDataset get_data() const;

            [[nodiscard]] std::vector<double> get_residuals(const std::vector<double>& params) override;

            /**
             * @brief Perform a fit and return the optimal parameters.
             */
            [[nodiscard]] std::vector<double> fit_params_only() override;

        protected:
            SimpleDataset data;
            std::unique_ptr<hist::DistanceHistogram> model;
            std::vector<mini::Parameter> guess;

            /**
			 * @brief Splice values from the model to match the data.
			 * 
			 * @param ym the model y-values corresponding to xm
			 */
			[[nodiscard]] std::vector<double> splice(const std::vector<double>& ym) const;

            /**
             * @brief Get the default guess parameters.
             */
            [[nodiscard]] std::vector<mini::Parameter> get_default_guess() const;

            /**
             * @brief Prepare a linear fitter for the given parameters.
             */
            detail::LinearLeastSquares prepare_linear_fitter(const std::vector<double>& params);
    };
}