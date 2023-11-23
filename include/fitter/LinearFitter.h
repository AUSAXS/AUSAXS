#pragma once

#include <fitter/Fitter.h>
#include <utility/observer_ptr.h>
#include <hist/HistFwd.h>
#include <io/IOFwd.h>

class Limit;
namespace fitter {

	/**
	 * @brief Fit an intensity curve to a dataset. 
	 * 
	 * Two parameters will be fitted: 
	 *    a: The slope of the curve.
	 *    b: The intercept of the curve.
	 * 
	 * This is just a convenient wrapper around SimpleLeastSquares. 
	 */
	class LinearFitter : public Fitter {
		protected: 
			LinearFitter() = default;

		public: 
            LinearFitter(LinearFitter&& other);

			/**
			 * Prepare a fit of the measured values in @a input to a model to be defined later. 
			 * 
			 * @param input The path to the file containing the measured values. 
			 */
			LinearFitter(const io::ExistingFile& input);

			/**
			 * Prepare a fit of the measured values in @a input to the model described by @a h.
			 * 
			 * @param input The path to the file containing the measured values. 
			 * @param h The distance histogram to fit. 
			 */
			LinearFitter(const io::ExistingFile& input, std::unique_ptr<hist::DistanceHistogram> h);

			/**
			 * Prepare a fit to the dataset.
			 */
			LinearFitter(const SimpleDataset& data);

			/**
			 * Prepare a fit of the histogram to the dataset. 
			 */
			LinearFitter(const SimpleDataset& data, std::unique_ptr<hist::DistanceHistogram> hist);

			/**
			 * Prepare a fit of the first histogram to the second. A series of data points is extracted from @a h2 and used to fit @a h1.
			 * 
			 * @param data The data histogram. 
			 * @param model The model histogram. 
			 */
			LinearFitter(std::unique_ptr<hist::DistanceHistogram> data, std::unique_ptr<hist::DistanceHistogram> model);

			/**
			 * Prepare a fit of the first histogram to the second. A series of data points is extracted from @a h2 and used to fit @a h1.
			 * 
			 * @param data The data histogram. 
			 * @param model The model histogram. 
			 * @param limits The limits on the generated data points. 
			 */
			LinearFitter(std::unique_ptr<hist::DistanceHistogram> data, std::unique_ptr<hist::DistanceHistogram> model, const Limit& limits);

			/**
			 * Prepare a fit to the histogram. A series of data points is extracted from it and used as the data points of the model. 
			 * 
			 * @param model The model histogram. 
			 */
			LinearFitter(std::unique_ptr<hist::DistanceHistogram> model);

			/**
			 * Prepare a fit to the histogram. A series of data points is extracted from it and used as the data points of the model. 
			 * 
			 * @param model The model histogram. 
			 * @param limits The limits on the generated data points. 
			 */
			LinearFitter(std::unique_ptr<hist::DistanceHistogram> model, const Limit& limits);

			/**
			 * @brief Destructor.
			 */
			virtual ~LinearFitter() override = default;

			/**
			 * @brief Perform the fit.
			 * 
			 * ! This function does NOT use the chi2 method, and is therefore not compatible with constraints.
			 * 
			 * @return A Fit object containing various information about the fit. Note that the fitted scaling parameter is a = c/M*r_e^2 and b = background
			 */
			[[nodiscard]] virtual std::shared_ptr<Fit> fit() override;

			/**
			 * @brief Perform the fit.
			 * 
			 * @return A Fit object containing various information about the fit. Note that the fitted scaling parameter is a = c/M*r_e^2 and b = background
			 */
            [[nodiscard]] virtual double fit_chi2_only() override;

			/**
			 * @brief Make a plot of the fit. 
			 * 
			 * @return A vector of TGraphs {Interpolated points, Optimal line, Measured points with uncertainties}
			 */
			[[nodiscard]] virtual FitPlots plot() override;

			/**
			 * @brief Make a residual plot of the fit.
			 * 
			 * @return A TGraphErrors with the residuals and their uncertainties. 
			 */
			[[nodiscard]] virtual SimpleDataset plot_residuals() override;

			/**
			 * @brief Change the scattering histogram used for the fit. 
			 */
			void set_scattering_hist(std::unique_ptr<hist::DistanceHistogram> h);

			/**
			 * @brief Get a view of the scattering histogram used for the fit. 
			 */
			[[nodiscard]] std::observer_ptr<hist::DistanceHistogram> get_scattering_hist();

			/**
			 * @brief Normalize all internally calculated intensities such that they start at this value.  
			 */
			void normalize_intensity(double I0);

			/**
			 * @brief Get the number of degrees of freedom. 
			 */
			[[nodiscard]] unsigned int degrees_of_freedom() const;

			/**
			 * @brief Get the number of degrees of freedom. 
			 */
			[[nodiscard]] unsigned int dof() const override;

			/**
			 * @brief Get the result of the last fit() call. 
			 */
			[[nodiscard]] virtual std::shared_ptr<Fit> get_fit() const override;

			void operator=(LinearFitter&& other);

		protected: 
			std::shared_ptr<Fit> fitted; 				// The previous fit result
			SimpleDataset data;          				// Observed data set
			double I0 = -1;              				// Normalization intensity
			std::unique_ptr<hist::DistanceHistogram> h; // The scattering histogram to fit

			/**
			 * @brief Calculate chi2 for a given choice of parameters @a params.
			 */
			[[nodiscard]] virtual double chi2(const std::vector<double>& params) override;

			/**
			 * @brief Prepare this class for fitting.
			 * 
			 * @param file measured values to compare the model against.
			 */
			void setup(const io::ExistingFile& file);

			/**
			 * @brief Splice values from the model to fit the evaluation points defined by the q values of the input file. 
			 * 
			 * @param ym the model y-values corresponding to xm
			 */
			[[nodiscard]] std::vector<double> splice(const std::vector<double>& ym) const;

			/**
			 * @brief Initialize this class based on a model histogram. 
			 */
			void model_setup(std::unique_ptr<hist::DistanceHistogram> model, const Limit& limits);
	};
}