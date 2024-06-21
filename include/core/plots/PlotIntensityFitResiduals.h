#pragma once

#include <plots/Plot.h>
#include <fitter/FitResult.h>
#include <utility/observer_ptr.h>

class SimpleDataset;
namespace plots {
	/**
	 * @brief Plot the residuals of the fitted scattering curve. 
	 * Remember to set the correct ScatteringPlot with the optimized values in the fitter before using this class. 
	 */
	class PlotIntensityFitResiduals : public Plot {
		public:
			/**
			 * @brief Constructor.
			 * 
			 * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
			 */
			PlotIntensityFitResiduals(fitter::LinearFitter& fitter);

			/**
			 * @brief Constructor.
			 * 
			 * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
			 */
			PlotIntensityFitResiduals(const fitter::FitResult& fit);

			/**
			 * @brief Constructor.
			 * 
			 * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
			 */
			PlotIntensityFitResiduals(observer_ptr<fitter::FitResult> fit);

			/**
			 * @brief Destructor.
			 */
			~PlotIntensityFitResiduals() override;

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(observer_ptr<fitter::FitResult> fit, const io::File& path);

		private:
			void plot(const SimpleDataset& graph);
	};
}