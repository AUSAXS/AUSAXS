#pragma once

#include <fitter/Fit.h>
#include <plots/Plot.h>
#include <histogram/ScatteringHistogram.h>
#include <fitter/SimpleIntensityFitter.h>

#include <memory>
#include <string>

namespace plots {

	/**
	 * @brief \class PlotIntensityFit
	 * 
	 * Plot both the measured and fitted scattering curve. 
	 * Remember to set the correct ScatteringPlot with the optimized values in the fitter before using this class. 
	 */
	class PlotIntensityFit : public Plot {
		public:
			/**
			 * @brief Constructor.
			 * 
			 * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
			 */
			PlotIntensityFit(SimpleIntensityFitter& fitter);

			/**
			 * @brief Constructor.
			 * 
			 * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
			 */
			PlotIntensityFit(const Fit& fit);

			/**
			 * @brief Constructor.
			 * 
			 * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
			 */
			PlotIntensityFit(const std::shared_ptr<Fit> fit);

			/**
			 * @brief Destructor.
			 */
			~PlotIntensityFit() override;

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(const std::shared_ptr<Fit> fit, std::string path);

			/**
			 * @brief Create and save the plot at the given path. 
			 * 
			 * @param path Save location and format. 
			 */
			void save(std::string path) const override;

		private:
			std::shared_ptr<TCanvas> canvas;

			void plot(const Fit::Plots& graphs) const;

			void prepare_canvas();
	};
}