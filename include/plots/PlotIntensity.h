#pragma once

#include <plots/Plot.h>
#include <plots/Styles.h>

#include <hist/HistFwd.h>
#include <fitter/FitterFwd.h>
#include <dataset/DatasetFwd.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace plots {
	class PlotIntensity : public Plot {
		public:
			PlotIntensity() = default;

			/**
			 * @brief Copy constructor.
			 * 
			 * @param d The ScatteringHistogram to be plotted. 
			 */
			PlotIntensity(const hist::ScatteringProfile& d, const plots::PlotOptions& options);

			/**
			 * @brief Constructor.
			 * 
			 * @param d The dataset to be plotted.
			 */
			PlotIntensity(const SimpleDataset& d, const plots::PlotOptions& options);

			/**
			 * @brief Destructor.
			 */
			~PlotIntensity() override;

			/**
			 * @brief Plot a scattering histogram.
			 */
			PlotIntensity& plot(const hist::ScatteringProfile& data, const plots::PlotOptions& options);

			/**
			 * @brief Plot an additional data set as points. 
			 */
			PlotIntensity& plot(const SimpleDataset& data, const plots::PlotOptions& options);

			/**
			 * @brief Plot the result of a fit. 
			 */
			PlotIntensity& plot(observer_ptr<fitter::Fit> fit, const plots::PlotOptions& options);

			/**
			 * @brief Plot the Guinier approximation for this scattering histogram. 
			 */
			PlotIntensity& plot_guinier_approx(observer_ptr<hist::ICompositeDistanceHistogram> data);

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(const hist::ScatteringProfile& h, const plots::PlotOptions& options, const io::File& path);

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(const SimpleDataset& d, const plots::PlotOptions& options, const io::File& path);
	};
}