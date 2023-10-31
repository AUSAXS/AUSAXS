#pragma once

#include <fitter/HydrationFitter.h>

#include <io/IOFwd.h>
#include <data/DataFwd.h>

namespace fitter {
    class Fit;
    class FitPlots;

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

            /**
             * @brief Destructor.
             */
            ~ExcludedVolumeFitter() override = default;

            /**
             * @brief Perform the fit.
             * 
             * @return A Fit object containing various information about the fit. Note that the fitted scaling parameter is a = c/M*r_e^2 and b = background
             */
            [[nodiscard]] std::shared_ptr<Fit> fit() override;

            [[nodiscard]] double fit_chi2_only() override;

            template<mini::type t>
            [[nodiscard]] std::shared_ptr<Fit> fit() {
                fit_type = t;
                return fit();
            }

            /**
             * @brief Make a plot of the fit. 
             * 
             * @return A vector of TGraphs {Interpolated points, Optimal line, Measured points with uncertainties}
             */
            [[nodiscard]] FitPlots plot() override;

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
             * @brief Set the guess value for the hydration scaling factor @a c.
             */
            void set_guess(std::vector<mini::Parameter> guess);

        private: 
            std::vector<mini::Parameter> guess; // The guess values for the parameters.
            mini::type fit_type;                // The algorithm to use.

            void update_excluded_volume(double d);

            /**
             * @brief Calculate chi2 for a given choice of parameters @a params.
             */
            [[nodiscard]] double chi2(const std::vector<double>& params) override;
    };
}