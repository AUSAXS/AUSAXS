#pragma once

#include <plots/Plot.h>
#include <fitter/FitResult.h>
#include <utility/observer_ptr.h>

namespace plots {
	/**
	 * @brief Plot both the measured and fitted scattering curve. 
	 * Remember to set the correct ScatteringPlot with the optimized values in the fitter before using this class. 
	 */
	class PlotIntensityFit : public Plot {
		public:
			/**
			 * @brief Constructor.
			 * 
			 * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
			 */
			PlotIntensityFit(fitter::LinearFitter& fitter);

			/**
			 * @brief Constructor.
			 * 
			 * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
			 */
			PlotIntensityFit(const fitter::FitResult& fit);

			/**
			 * @brief Constructor.
			 * 
			 * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
			 */
			PlotIntensityFit(observer_ptr<fitter::FitResult> fit);

			/**
			 * @brief Destructor.
			 */
			~PlotIntensityFit() override;

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(observer_ptr<fitter::FitResult> fit, const io::File& path);

		private:
			void plot(const fitter::FitResult::FitPlots& graphs);
	};
}