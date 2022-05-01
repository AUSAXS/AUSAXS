#pragma once

#include <vector>
#include <string>
#include <memory>

#include <fitter/Fit.h>
#include <fitter/Fitter.h>

/**
 * @brief A simple linear least-squares fitter for fitting the linear relationship y = ax+b.
 */
class SimpleLeastSquares : public Fitter {
  public:
    /**
     * @brief Constructor.
     * 
     * Prepare a linear least-squares fit for the given dataset. 
     */
    SimpleLeastSquares(const Dataset& data);

    /**
     * @brief Constructor.
     * 
     * Prepare a linear least-squares fit for the inputs. 
     */
    SimpleLeastSquares(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& y_err);

    /**
     * @brief Destructor.
     */
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

    Multiset plot() override;

    Dataset plot_residuals() override;

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