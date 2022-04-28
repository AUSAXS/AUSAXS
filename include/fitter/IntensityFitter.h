#pragma once

#include <fitter/Fit.h>
#include <fitter/SimpleIntensityFitter.h>
#include <ScatteringHistogram.h>

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
     *        Prepare a fit of the histogram to the measured values. 
     * 
     * @param input The path to the file containing the measured values. 
     * @param h The histogram.
     */
    IntensityFitter(std::string input, const ScatteringHistogram& h) : SimpleIntensityFitter(input, h) {}

    /**
     * @brief Constructor.
     *        Prepare a fit of the histogram to the measured values. 
     * 
     * @param input The path to the file containing the measured values. 
     * @param h The histogram.
     */
    IntensityFitter(std::string input, ScatteringHistogram&& h) : SimpleIntensityFitter(input, h) {}

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

    /**
     * @brief Make a plot of the fit. 
     * 
     * @return A vector of TGraphs {Interpolated points, Optimal line, Measured points with uncertainties}
     */
    std::vector<std::shared_ptr<TGraph>> plot() override;

    /**
     * @brief Make a residual plot of the fit.
     * 
     * @return A TGraphErrors with the residuals and their uncertainties. 
     */
    std::unique_ptr<TGraphErrors> plot_residuals() override;

    /**
     * @brief Get the intercept of the model. This might be useful for calculating the concentration.
     */
    double get_intercept();

    /**
     * @brief Get the model dataset for the points specified by the input data file. 
     */
    SAXSDataset get_model_dataset();

    /**
     * @brief Get the model dataset for the points specified by @a q. 
     */
    SAXSDataset get_model_dataset(const vector<double>& q);

    /**
     * @brief Get the dataset being fitted. 
     */
    SAXSDataset get_dataset() const;

  private: 
    /**
     * @brief Calculate chi2 for a given choice of parameters @a params.
     */
    double chi2(const double* params) override;
};