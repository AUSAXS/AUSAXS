#pragma once

#include <fitter/Fit.h>
#include <fitter/Fitter.h>
#include <hist/ScatteringHistogram.h>

/**
 * @brief Perform a simple chi2 fit of a data set to a scattering curve. 
 * 
 * Two parameters will be fitted: 
 *    a: The slope of the curve.
 *    b: The intercept of the curve.
 * 
 * This is just a convenient wrapper around SimpleLeastSquares. 
 */
class SimpleIntensityFitter : public Fitter {
	public: 
		/**
		 * @brief Constructor.
		 * 
		 * Prepare a fit of the measured values in @a input to a model to be defined later. 
		 * 
		 * @param input The path to the file containing the measured values. 
		 */
		SimpleIntensityFitter(std::string input) {setup(input);}

		/**
		 * @brief Constructor.
		 * 
		 * Prepare a fit of the measured values in @a input to the model described by @a h.
		 * 
		 * @param input The path to the file containing the measured values. 
		 * @param h The ScatteringHistogram to fit. 
		 */
		SimpleIntensityFitter(std::string input, const hist::ScatteringHistogram& h) : h(h) {setup(input);}

		/**
		 * @brief Constructor.
		 * 
		 * Prepare a fit of the measured values in @a input to the model described by @a h.
		 * 
		 * @param input the path to the file containing the measured values. 
		 * @param h The ScatteringHistogram to fit. 
		 */
		SimpleIntensityFitter(std::string input, hist::ScatteringHistogram&& h) : h(std::move(h)) {setup(input);}

		/**
		 * @brief Constructor. 
		 * 
		 * Prepare a fit to the dataset.
		 */
		SimpleIntensityFitter(const SimpleDataset& data) : data(data) {}

		/**
		 * @brief Constructor. 
		 * 
		 * Prepare a fit of the histogram to the dataset. 
		 */
		SimpleIntensityFitter(const SimpleDataset& data, const hist::ScatteringHistogram& hist) : data(data), h(hist) {}

		/**
		 * @brief Constructor.
		 * 
		 * Prepare a fit of the first histogram to the second. A series of data points is extracted from @a h2 and used to fit @a h1.
		 * 
		 * @param data The data histogram. 
		 * @param model The model histogram. 
		 * @param limits The limits on the generated data points. 
		 */
		SimpleIntensityFitter(const hist::ScatteringHistogram& data, const hist::ScatteringHistogram& model, const Limit& limits = Limit(setting::axes::qmin, setting::axes::qmax));

		/**
		 * @brief Constructor.
		 * 
		 * Prepare a fit to the histogram. A series of data points is extracted from it and used as the data points of the model. 
		 * 
		 * @param model The model histogram. 
		 * @param limits The limits on the generated data points. 
		 */
		SimpleIntensityFitter(const hist::ScatteringHistogram& model, const Limit& limits = Limit(setting::axes::qmin, setting::axes::qmax));

		/**
		 * @brief Destructor.
		 */
		virtual ~SimpleIntensityFitter() override = default;

		/**
		 * @brief Perform the fit.
		 * 
		 * @return A Fit object containing various information about the fit. Note that the fitted scaling parameter is a = c/M*r_e^2 and b = background
		 */
		virtual std::shared_ptr<Fit> fit() override;

		/**
		 * @brief Make a plot of the fit. 
		 * 
		 * @return A vector of TGraphs {Interpolated points, Optimal line, Measured points with uncertainties}
		 */
		virtual Fit::Plots plot() override;

		/**
		 * @brief Make a residual plot of the fit.
		 * 
		 * @return A TGraphErrors with the residuals and their uncertainties. 
		 */
		virtual SimpleDataset plot_residuals() override;

		/**
		 * @brief Change the scattering histogram used for the fit. 
		 */
		void set_scattering_hist(const hist::ScatteringHistogram& h);

		/**
		 * @brief Change the scattering histogram used for the fit. 
		 */
		void set_scattering_hist(hist::ScatteringHistogram&& h);

		/**
		 * @brief Normalize all internally calculated intensities such that they start at this value.  
		 */
		void normalize_intensity(double I0);

		/**
		 * @brief Get the number of degrees of freedom. 
		 */
		unsigned int degrees_of_freedom() const;

		/**
		 * @brief Get the number of degrees of freedom. 
		 */
		unsigned int dof() const override;

		/**
		 * @brief Get the result of the last fit() call. 
		 */
		virtual std::shared_ptr<Fit> get_fit() const override;

	protected: 
		std::shared_ptr<Fit> fitted; // The previous fit result
		SimpleDataset data;          // Observed data set
		double I0 = -1;              // Normalization intensity
		hist::ScatteringHistogram h; // The scattering histogram to fit

		/**
		 * @brief Calculate chi2 for a given choice of parameters @a params.
		 */
		virtual double chi2(const double* params);

		/**
		 * @brief Prepare this class for fitting.
		 * 
		 * @param file measured values to compare the model against.
		 */
		void setup(std::string file);

		/**
		 * @brief Splice values from the model to fit the evaluation points defined by the q values of the input file. 
		 * 
		 * @param ym the model y-values corresponding to xm
		 */
		std::vector<double> splice(const std::vector<double>& ym) const;

		/**
		 * @brief Initialize this class based on a model histogram. 
		 */
		void model_setup(const hist::ScatteringHistogram& model, const Limit& limits);
};