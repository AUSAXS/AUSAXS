#pragma once

#include <fitter/LinearFitter.h>
#include <mini/detail/Parameter.h>
#include <io/IOFwd.h>

namespace fitter {
    class FitPlots;

    /**
     * @brief Fit an intensity curve to a dataset. 
     * 
     * Three parameters will be fitted: 
     *    a: The slope of the curve.
     *    b: The intercept of the curve.
     *    c: The scattering length of the hydration shell.
     */
    class HydrationFitter : public LinearFitter {
        protected:
            HydrationFitter() = default;

        public: 
            HydrationFitter(HydrationFitter&& other);

            /**
             * Prepare a fit of the measured values in @a input to a model to be defined later. 
             * 
             * @param input The path to the file containing the measured values. 
             */
            HydrationFitter(const io::ExistingFile& input);

            /**
             * Prepare a fit of the histogram to the measured values. 
             * 
             * @param input The path to the file containing the measured values. 
             * @param h The histogram.
             */
            HydrationFitter(const io::ExistingFile& input, std::unique_ptr<hist::ICompositeDistanceHistogram> h);

            /**
             * Prepare a fit of the histogram to the measured data. 
             * 
             * @param data The measured data.
             * @param h The histogram.
             */
            HydrationFitter(const SimpleDataset& data, std::unique_ptr<hist::ICompositeDistanceHistogram> h);

            /**
             * @brief Constructor.
             * 
             * Prepare a fit to the histogram. A series of data points is extracted from it and used as the data points of the model. 
             * 
             * @param model The model histogram. 
             * @param limits The limits on the generated data points. 
             */
            HydrationFitter(std::unique_ptr<hist::ICompositeDistanceHistogram> model);

            /**
             * @brief Constructor.
             * 
             * Prepare a fit to the histogram. A series of data points is extracted from it and used as the data points of the model. 
             * 
             * @param model The model histogram. 
             * @param limits The limits on the generated data points. 
             */
            HydrationFitter(std::unique_ptr<hist::ICompositeDistanceHistogram> model, const Limit& limits);

            /**
             * @brief Destructor.
             */
            virtual ~HydrationFitter() override = default;

            /**
             * @brief Perform the fit.
             * 
             * @return A Fit object containing various information about the fit. Note that the fitted scaling parameter is a = c/M*r_e^2 and b = background
             */
            [[nodiscard]] virtual std::shared_ptr<Fit> fit() override;

            [[nodiscard]] virtual double fit_chi2_only() override;

            [[nodiscard]] std::shared_ptr<Fit> fit(const mini::type& algorithm);

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
             *        These points are interpolated from the model histogram.
             */
            [[nodiscard]] SimpleDataset get_model_dataset(const std::vector<double>& q);

            /**
             * @brief Get the dataset being fitted. 
             */
            [[nodiscard]] SimpleDataset get_dataset() const;

            /**
             * @brief Set the guess value for the hydration scaling factor @a c.
             */
            void set_guess(const mini::Parameter& guess);

            /**
             * @brief Set the fitting algorithm to use.
             */
            void set_algorithm(const mini::type& t);

            void operator=(HydrationFitter&& other);

            /**
			 * @brief Get a view of the scattering histogram used for the fit. 
			 */
			[[nodiscard]] view_ptr<hist::ICompositeDistanceHistogram> get_scattering_hist();

            mini::type fit_type = mini::type::BFGS;
        protected:
            /**
             * @brief Calculate chi2 for a given choice of parameters @a params.
             */
            [[nodiscard]] virtual double chi2(const std::vector<double>& params) override;

            /**
             * @brief Cast the histogram to a CompositeDistanceHistogram.
             * 
             * @return This is always safe since the constructor of this class guarantees that the histogram is a CompositeDistanceHistogram.
             */
            view_ptr<hist::ICompositeDistanceHistogram> cast_h() const;

        private: 
            mini::Parameter guess = {"c", 1, {0, 10}}; // The guess value for the hydration scaling factor.
    };
}