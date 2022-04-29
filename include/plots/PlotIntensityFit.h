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
		 * @brief Create and save the plot at the given path. 
		 * 
		 * @param path Save location and format. 
		 */
		void save(std::string path) const override;

		private:
		std::unique_ptr<TCanvas> canvas;

		void plot(const std::vector<std::shared_ptr<TGraph>>& graphs) const;

		void prepare_canvas();
	};
}