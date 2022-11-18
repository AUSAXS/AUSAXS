#pragma once

#include <fitter/Fit.h>
#include <fitter/SimpleIntensityFitter.h>
#include <hist/ScatteringHistogram.h>

/**
 * @brief Fit an intensity curve to a dataset. 
 * 
 * Three parameters will be fitted: 
 *    a: The slope of the curve
 *    b: The intercept of the curve
 *    c: The scaling factor for the hydration shell
 */
class IntensityFitter : public SimpleIntensityFitter {
    public: 
        /**
         * @brief Constructor.
         * 
         * Prepare a fit of the measured values in @a input to a model to be defined later. 
         * 
         * @param input The path to the file containing the measured values. 
         */
        IntensityFitter(std::string input) : SimpleIntensityFitter(input) {}

        /**
         * @brief Constructor.
         *        Prepare a fit of the histogram to the measured values. 
         * 
         * @param input The path to the file containing the measured values. 
         * @param h The histogram.
         */
        IntensityFitter(std::string input, const hist::ScatteringHistogram& h) : SimpleIntensityFitter(input, h) {}

        /**
         * @brief Constructor.
         *        Prepare a fit of the histogram to the measured values. 
         * 
         * @param input The path to the file containing the measured values. 
         * @param h The histogram.
         */
        IntensityFitter(std::string input, hist::ScatteringHistogram&& h) : SimpleIntensityFitter(input, h) {}

        /**
         * @brief Constructor.
         *        Prepare a fit of the histogram to the measured data. 
         * 
         * @param data The measured data.
         * @param h The histogram.
         */
		IntensityFitter(const SimpleDataset& data, const hist::ScatteringHistogram& h) : SimpleIntensityFitter(data, h) {}

        /**
         * @brief Constructor.
         * 
         * Prepare a fit to the histogram. A series of data points is extracted from it and used as the data points of the model. 
         * 
         * @param model The model histogram. 
         * @param limits The limits on the generated data points. 
         */
        IntensityFitter(const hist::ScatteringHistogram& model, const Limit& limits = Limit(setting::axes::qmin, setting::axes::qmax));

        /**
         * @brief Destructor.
         */
        ~IntensityFitter() override = default;

        /**
         * @brief Perform the fit.
         * 
         * @return A Fit object containing various information about the fit. Note that the fitted scaling parameter is a = c/M*r_e^2 and b = background
         */
        std::shared_ptr<Fit> fit() override;

        template<mini::Type t>
        std::shared_ptr<Fit> fit() {
            fit_type = t;
            return fit();
        }

        /**
         * @brief Make a plot of the fit. 
         * 
         * @return A vector of TGraphs {Interpolated points, Optimal line, Measured points with uncertainties}
         */
        Fit::Plots plot() override;

        /**
         * @brief Make a residual plot of the fit.
         * 
         * @return A TGraphErrors with the residuals and their uncertainties. 
         */
        SimpleDataset plot_residuals() override;

        /**
         * @brief Get the intercept of the model. This might be useful for calculating the concentration.
         */
        double get_intercept();

        /**
         * @brief Get the model dataset for the points specified by the input data file. 
         */
        SimpleDataset get_model_dataset();

        /**
         * @brief Get the model dataset for the points specified by @a q. 
         */
        SimpleDataset get_model_dataset(const std::vector<double>& q);

        /**
         * @brief Get the dataset being fitted. 
         */
        SimpleDataset get_dataset() const;

        /**
         * @brief Set the guess value for the hydration scaling factor @a c.
         */
        void set_guess(mini::Parameter guess);

    private: 
    	mini::Parameter guess = {"c", 5, {0, 10}}; // The guess value for the hydration scaling factor.
        mini::Type fit_type = mini::Type::BFGS;    // The algorithm to use.

        /**
         * @brief Calculate chi2 for a given choice of parameters @a params.
         */
        double chi2(std::vector<double> params) override;
};