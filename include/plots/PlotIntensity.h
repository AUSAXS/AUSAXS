#pragma once

#include <plots/Plot.h>
#include <plots/Styles.h>

#include <memory>

namespace hist {class ScatteringHistogram;}
namespace fitter {
	class Fitter;
	class Fit;
}
class SimpleDataset;
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
			PlotIntensity& plot(const hist::ScatteringHistogram& data, style::Color color = style::color::black);

			/**
			 * @brief Plot an additional data set as points. 
			 */
			PlotIntensity& plot(const SimpleDataset& data, style::Color color = style::color::black);

			/**
			 * @brief Plot the result of a fit. 
			 */
			PlotIntensity& plot(const std::shared_ptr<fitter::Fit> fit, style::Color color = style::color::black);

			/**
			 * @brief Plot the Guinier approximation for this scattering histogram. 
			 */
			PlotIntensity& plot_guinier_approx(const hist::ScatteringHistogram& data);

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(const hist::ScatteringHistogram& h, const io::File& path);

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(const SimpleDataset& d, const io::File& path);
	};
}