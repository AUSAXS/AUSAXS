#pragma once

class TPad;

#include <fitter/Fit.h>
#include <fitter/Fitter.h>
#include <plots/Plot.h>
#include <histogram/ScatteringHistogram.h>

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
		PlotIntensity(const ScatteringHistogram& d, int color = kBlack);

		/**
		 * @brief Move constructor.
		 * 
		 * @param d The ScatteringHistogram to be plotted. 
		 */
		PlotIntensity(ScatteringHistogram&& d, int color = kBlack);

		/**
		 * @brief Constructor.
		 * 
		 * @param d The dataset to be plotted.
		 */
		PlotIntensity(const SAXSDataset& d);

		/**
		 * @brief Destructor.
		 */
		~PlotIntensity() override;

		/**
		 * @brief Plot an additional data set as points. 
		 */
		void plot_intensity(const Dataset& data);

		/**
		 * @brief Plot the result of a fit. 
		 */
		void plot_intensity(const std::shared_ptr<Fit> fit, const PlotOptions& options);

		void plot_guinier_approx();

		void save(std::string path) const override;

		private:
		const ScatteringHistogram d;
		std::unique_ptr<TCanvas> canvas;
		std::unique_ptr<TPad> linpad;
		std::unique_ptr<TPad> logpad;
		double ymin, ymax;

		void initial_intensity_plot(int color);

		void initial_intensity_plot(const Dataset& data);

		void prepare_canvas();
	};
}