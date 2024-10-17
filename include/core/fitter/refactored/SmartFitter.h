#pragma once

#include <fitter/refactored/Fitter.h>
#include <dataset/SimpleDataset.h>
#include <mini/detail/Parameter.h>
#include <hist/HistFwd.h>

namespace fitter {
    /**
     * @brief A smart fitter automatically fitting the parameters enabled in the settings. 
     *
     * Supported parameters: linear (always), solvent shell density, excluded volume, solvent density
     */
    class SmartFitter : public Fitter {
        public:
            /**
             * @brief Prepare a fit of the measured values in @a input to a model to be defined later. 
             */
            SmartFitter(const SimpleDataset& data);

            /**
             * @brief Prepare a fit of the measured values in @a input to the model described by @a h.
             */
            SmartFitter(const SimpleDataset& saxs, std::unique_ptr<hist::ICompositeDistanceHistogram> h);
            virtual ~SmartFitter() override;

            [[nodiscard]] virtual std::unique_ptr<FitResult> fit() override;
            [[nodiscard]] virtual double fit_chi2_only() override;
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
             * @brief Get the optimized model. 
             *
             * @throw bad_order if the fit has not been performed yet.
             */
            SimpleDataset get_optimized_model() const;

            /**
             * @brief Get the dataset being fitted. 
             */
            SimpleDataset get_data() const;

        private:
            SimpleDataset data;
            std::unique_ptr<hist::ICompositeDistanceHistogram> model;
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

            [[nodiscard]] double chi2(const std::vector<double>& params) override;
    };
}