#pragma once

#include <fitter/Fit.h>
#include <fitter/Fitter.h>
#include <plots/Plot.h>
#include <hist/ScatteringHistogram.h>

#include <memory>
#include <string>

namespace plots {
	class PlotIntensity : public Plot {
		public:
			/**
			 * @brief Copy constructor.
			 * 
			 * @param d The ScatteringHistogram to be plotted. 
			 */
			PlotIntensity(const hist::ScatteringHistogram& d, style::Color color = style::color::black);

			/**
			 * @brief Constructor.
			 * 
			 * @param d The dataset to be plotted.
			 */
			PlotIntensity(const SimpleDataset& d, style::Color color = style::color::black);

			/**
			 * @brief Destructor.
			 */
			~PlotIntensity() override;

			/**
			 * @brief Plot a scattering histogram.
			 */
			void plot(const hist::ScatteringHistogram& data, style::Color color = style::color::black);

			/**
			 * @brief Plot an additional data set as points. 
			 */
			void plot(const SimpleDataset& data, style::Color color = style::color::black);

			/**
			 * @brief Plot the result of a fit. 
			 */
			void plot(const std::shared_ptr<Fit> fit, style::Color color = style::color::black);

			/**
			 * @brief Plot the Guinier approximation for this scattering histogram. 
			 */
			void plot_guinier_approx(const hist::ScatteringHistogram& data);

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(const hist::ScatteringHistogram& h, std::string path);

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(const SimpleDataset& d, std::string path);
	};
}