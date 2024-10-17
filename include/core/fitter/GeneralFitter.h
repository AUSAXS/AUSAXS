#pragma once

#include <fitter/HydrationFitter.h>

#include <io/IOFwd.h>
#include <data/DataFwd.h>
#include <fitter/FitterFwd.h>

namespace fitter {
    class AutoFitter : public LinearFitter {
        public:
            AutoFitter(const SimpleDataset& data);
            AutoFitter(const SimpleDataset& saxs, std::unique_ptr<hist::ICompositeDistanceHistogram> h);
            ~AutoFitter() override = default;

            /**
             * @brief Perform the fit.
             * 
             * @return A Fit object containing various information about the fit. 
             */
            [[nodiscard]] std::shared_ptr<FitResult> fit() override;
            [[nodiscard]] double fit_chi2_only() override;

            [[nodiscard]] unsigned int dof() const override;

            /**
             * @brief Set the guess values for the fit. 
             */
            void set_guess(std::vector<mini::Parameter> guess);

            /**
             * @brief Set the scattering histogram to use for the fit. 
             */
            void set_scattering_hist(std::unique_ptr<hist::ICompositeDistanceHistogram> h);

        private: 
            /**
             * @brief Validate that the histogram is compatible with the current fitter.
             * 
             * @throw invalid_argument if the histogram is not compatible.
             */
            void validate_histogram() const;

            /**
             * @brief Calculate chi2 for a given choice of parameters @a params.
             */
            [[nodiscard]] double chi2(const std::vector<double>& params) override;
    };
}