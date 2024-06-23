#pragma once

#include <fitter/HydrationFitter.h>

#include <io/IOFwd.h>
#include <data/DataFwd.h>
#include <fitter/FitterFwd.h>

namespace fitter {
    /**
     * @brief Fit an intensity curve to a dataset. 
     * 
     * Four parameters will be fitted: 
     *    a: The slope of the curve.
     *    b: The intercept of the curve.
     *    c: The scattering length of the hydration shell.
     *    d: The excluded volume. 
     */
    class ExcludedVolumeFitter : public HydrationFitter {
        public:
            ExcludedVolumeFitter(const io::ExistingFile& saxs, std::unique_ptr<hist::ICompositeDistanceHistogram> h);

            ~ExcludedVolumeFitter() override = default;

            /**
             * @brief Perform the fit.
             * 
             * @return A Fit object containing various information about the fit. Note that the fitted scaling parameter is a = c/M*r_e^2 and b = background
             */
            [[nodiscard]] std::shared_ptr<FitResult> fit() override;

            [[nodiscard]] double fit_chi2_only() override;

            template<mini::type t>
            [[nodiscard]] std::shared_ptr<FitResult> fit() {
                fit_type = t;
                return fit();
            }

            /**
             * @brief Make a plot of the fit. 
             * 
             * @return A vector of TGraphs {Interpolated points, Optimal line, Measured points with uncertainties}
             */
            [[nodiscard]] FitResult::FitInfo plot() override;

            /**
             * @brief Make a residual plot of the fit.
             * 
             * @return A TGraphErrors with the residuals and their uncertainties. 
             */
            [[nodiscard]] SimpleDataset plot_residuals() override;

            /**
             * @brief Get the intercept of the model. This might be useful for calculating the concentration.
             */
            [[nodiscard]] double get_intercept();

            /**
             * @brief Get the model dataset for the points specified by the input data file. 
             */
            [[nodiscard]] SimpleDataset get_model_dataset();

            /**
             * @brief Get the model dataset for the points specified by @a q. 
             */
            [[nodiscard]] SimpleDataset get_model_dataset(const std::vector<double>& q);

            /**
             * @brief Get the dataset being fitted. 
             */
            [[nodiscard]] SimpleDataset get_dataset() const;

            /**
             * @brief Set the guess values for the hydration scaling factor @a c.
             */
            void set_guess(mini::Parameter guess_hydration, mini::Parameter guess_exv);

            /**
             * @brief Set the scattering histogram to use for the fit. 
             */
            void set_scattering_hist(std::unique_ptr<hist::ICompositeDistanceHistogram> h);

        private: 
            void initialize_guess() override;
            mini::Parameter guess_exv;  // The guess values for the parameters.

            void update_excluded_volume(double d);

            /**
             * @brief Validate that the histogram is a ICompositeDistanceHistogramExv.
             *        Since the ICompositeDistanceHistogramExv class is rarely used elsewhere, it cannot be used directly in the constructor.
             */
            void validate_histogram() const;

            /**
             * @brief Cast the histogram to a ICompositeDistanceHistogramExv.
             * 
             * @return This is always safe since the constructor of this class guarantees that the histogram is a ICompositeDistanceHistogramExv.
             */
            observer_ptr<hist::ICompositeDistanceHistogramExv> cast_exv() const;

            /**
             * @brief Calculate chi2 for a given choice of parameters @a params.
             */
            [[nodiscard]] double chi2(const std::vector<double>& params) override;
    };
}