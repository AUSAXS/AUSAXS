#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <memory>

#include "fitter/Fitter.h"

#include <Math/SpecFuncMathCore.h> // for the incomplete gamma function

using std::vector, std::shared_ptr;

/**
 * @brief A simple linear least-squares fitter for fitting the linear relationship y = ax+b.
 */
class SimpleLeastSquares : public Fitter {
    public:
        SimpleLeastSquares(const vector<double>& x, const vector<double>& y, const vector<double>& y_err) : x(x), y(y), y_err(y_err) {}
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
        virtual shared_ptr<Fit> fit() override;

    private:
        const vector<double> &x, &y, &y_err;
        double S, Sx, Sy, Sxx, Sxy, delta = 0;
        double a, b;

        /**
         * @brief Calculate chi2.
         */
        double chi2() const;
};