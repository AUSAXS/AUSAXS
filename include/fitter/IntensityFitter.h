#pragma once

#include <fitter/Fit.h>
#include <fitter/SimpleIntensityFitter.h>
#include <ScatteringHistogram.h>

class IntensityFitter : public SimpleIntensityFitter {
  public: 
    /**
     * @brief Constructor.
     *        Prepare a fit of the measured values in @a input to the model described by @a q and @a I.
     * 
     * @param input the path to the file containing the measured values. 
     * @param q the model q values.
     * @param I the model I values. 
     */
    IntensityFitter(std::string input, const ScatteringHistogram& h) : SimpleIntensityFitter(input, h) {}

    /**
     * @brief Constructor.
     *        Prepare a fit of the measured values in @a input to the model described by @a q and @a I.
     * 
     * @param input the path to the file containing the measured values. 
     * @param q the model q values.
     * @param I the model I values. 
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

  private: 
    /**
     * @brief Calculate chi2 for a given choice of parameters @a params.
     */
    double chi2(const double* params) override;
};