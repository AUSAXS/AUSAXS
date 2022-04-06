#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <memory>

#include "fitter/Fitter.h"

#include <Math/SpecFuncMathCore.h> // for the incomplete gamma function

/**
 * @brief A simple linear least-squares fitter for fitting the linear relationship y = ax+b.
 */
class SimpleLeastSquares : public Fitter {
  public:
    SimpleLeastSquares(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& y_err) : x(x), y(y), y_err(y_err) {}
    ~SimpleLeastSquares() override {}

    /**
     * @brief Perform a linear least-squares fit and calculate @a only the fitted parameters.
     *        There are no guarantees on the goodness-of-fit. Use the standard @a fit() instead for that. 
     * @return The fitted parameters (a, b) for the equation y = ax+b.
     */
    std::pair<double, double> fit_params_only();

    /**
     * @brief Perform a linear least-squares fit. 
     * @return A Fit object containing various information for the fit. 
     */
    virtual std::shared_ptr<Fit> fit() override;

    std::vector<std::shared_ptr<TGraph>> plot() override;

    std::unique_ptr<TGraphErrors> plot_residuals() override;

    unsigned int dof() const override;

  private:
    const std::vector<double> &x, &y, &y_err;
    double S, Sx, Sy, Sxx, Sxy, delta = 0;
    double a, b;

    /**
     * @brief Calculate chi2.
     */
    double chi2() const;
};