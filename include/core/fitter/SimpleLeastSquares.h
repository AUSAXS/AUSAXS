#pragma once

#include <fitter/Fit.h>
#include <fitter/Fitter.h>

#include <vector>
#include <memory>

namespace fitter {
    /**
     * @brief A simple linear least-squares fitter for fitting the linear relationship y = ax+b.
     */
    class SimpleLeastSquares : public Fitter {
        public:
            /**
             * @brief Prepare a linear least-squares fit for the given dataset. 
             */
            SimpleLeastSquares(const std::vector<double> data, const std::vector<double> model);

            /**
             * @brief Prepare a linear least-squares fit for the given dataset. 
             */
            SimpleLeastSquares(const std::vector<double> data, const std::vector<double> model, const std::vector<double> errors);

            /**
             * @brief Prepare a linear least-squares fit for the given dataset. 
             */
            SimpleLeastSquares(const SimpleDataset& data);

            /**
             * @brief Prepare a linear least-squares fit for the given dataset. 
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
            [[nodiscard]] std::pair<double, double> fit_params_only();

            [[nodiscard]] virtual double fit_chi2_only() override;

            /**
             * @brief Perform a linear least-squares fit. 
             * @return A Fit object containing various information for the fit. 
             */
            [[nodiscard]] virtual std::shared_ptr<Fit> fit() override;

            /**
             * @brief Get a multiset containing the fitted curve of the last fit() call. 
             */
            [[nodiscard]] FitPlots plot() override;

            /**
             * @brief Get a dataset containing the residuals of the last fit() call. 
             */
            [[nodiscard]] SimpleDataset plot_residuals() override;

            /**
             * @brief Get the result of the last fit() call.
             */
            [[nodiscard]] virtual std::shared_ptr<Fit> get_fit() const override; 

            /**
             * @brief Get the number of degrees of freedom.
             */
            [[nodiscard]] unsigned int dof() const override;

            [[nodiscard]] unsigned int size() const override;

        private:
            const SimpleDataset data;
            double S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0, delta = 0;
            double a = 0, b = 0;

            /**
             * @brief Calculate chi2.
             */
            [[nodiscard]] double chi2(const std::vector<double>& params) override;
    };
}