#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <memory>

#include "Fitter.h"

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
        std::pair<double, double> fit_params_only() {
            S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
            for (size_t i = 0; i < x.size(); i++) {
                double sig2 = pow(y_err[i], 2);
                S += 1./sig2;
                Sx += x[i]/sig2;
                Sy += y[i]/sig2;
                Sxx += pow(x[i], 2)/sig2;
                Sxy += x[i]*y[i]/sig2;
            }

            delta = S*Sxx - pow(Sx, 2);
            a = (S*Sxy - Sx*Sy)/delta;
            b = (Sxx*Sy - Sx*Sxy)/delta;
            return std::make_pair(a, b);
        }

        /**
         * @brief Perform a linear least-squares fit. 
         * @return A Fit object containing various information for the fit. 
         */
        virtual shared_ptr<Fit> fit() override {
            if (delta == 0) {fit_params_only();}
            double a_err2 = S/delta; // squared sigmas
            double b_err2 = Sxx/delta; 

            // double cov_ab = -Sx/delta;
            double Q = ROOT::Math::inc_gamma((double) x.size()/2 -1, chi2()/2);

            shared_ptr<Fit> f = std::make_shared<Fit>();
            f->params = {{"a", a}, {"b", b}};
            f->errs = {{"a", sqrt(a_err2)}, {"b", sqrt(b_err2)}};
            f->dof = x.size() - 2;
            f->chi2 = chi2();
            f->calls = 1;
            f->converged = Q > 0.001;
            return f;
        }

    private:
        const vector<double> &x, &y, &y_err;
        double S, Sx, Sy, Sxx, Sxy, delta = 0;
        double a, b;

        /**
         * @brief Calculate chi2.
         */
        double chi2() const {
            double chi = 0;
            for (size_t i = 0; i < x.size(); ++i) {
                chi += pow((y[i] - a*x[i] - b)/y_err[i], 2);
            }
            return chi;
        }
};