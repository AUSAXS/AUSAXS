#pragma once

#include <plots/Plot.h>

#include <memory>

namespace fitter {
	class LinearFitter;
	class Fit;
}
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
			PlotIntensityFitResiduals(const fitter::Fit& fit);

			/**
			 * @brief Constructor.
			 * 
			 * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
			 */
			PlotIntensityFitResiduals(const std::shared_ptr<fitter::Fit> fit);

			/**
			 * @brief Destructor.
			 */
			~PlotIntensityFitResiduals() override;

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single histogram. 
			 */
			static void quick_plot(const std::shared_ptr<fitter::Fit> fit, const io::File& path);

		private:
			void plot(const SimpleDataset& graph);
	};
}