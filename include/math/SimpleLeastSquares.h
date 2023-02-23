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
         * @brief Copy constructor.
         * 
         * Prepare a linear least-squares fit for the given dataset. 
         */
        SimpleLeastSquares(const SimpleDataset& data);

        /**
         * @brief Move constructor.
         * 
         * Prepare a linear least-squares fit for the given dataset. 
         */
        SimpleLeastSquares(SimpleDataset&& data);

        /**
         * @brief Destructor.
         */
        virtual ~SimpleLeastSquares() override = default;

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

        /**
         * @brief Get a multiset containing the fitted curve of the last fit() call. 
         */
        Fit::Plots plot() override;

        /**
         * @brief Get a dataset containing the residuals of the last fit() call. 
         */
        SimpleDataset plot_residuals() override;

        /**
         * @brief Get the result of the last fit() call.
         */
        virtual std::shared_ptr<Fit> get_fit() const override; 

        /**
         * @brief Get the number of degrees of freedom.
         */
        unsigned int dof() const override;

    private:
        const SimpleDataset data;
        double S, Sx, Sy, Sxx, Sxy, delta = 0;
        double a, b;

        /**
         * @brief Calculate chi2.
         */
        [[nodiscard]] double chi2(const std::vector<double>& params) override;
};