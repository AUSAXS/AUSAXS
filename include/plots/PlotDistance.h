#pragma once

#include <plots/Plot.h>
#include <histogram/ScatteringHistogram.h>

namespace plots {
	/**
	 * @brief \class PlotDistance.
	 * 
	 * Plots a histogram of all distances. 
	 */
	class PlotDistance : public Plot {
		public:
			/**
			 * @brief Constructor.
			 * 
			 * @param d The ScatteringHistogram which will be plotted. 
			 * @param path The path to the folder where the plot will be saved. 
			 */
			PlotDistance(const hist::ScatteringHistogram& d, std::string path);

			/**
			 * @brief Destructor. 
			 */
			~PlotDistance() override;
		
		private: 
			void plot(const hist::ScatteringHistogram& d);
	};
}