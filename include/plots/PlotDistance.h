#pragma once

#include <plots/Plot.h>
#include <hist/ScatteringHistogram.h>

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
			 * @param data The ScatteringHistogram which will be plotted. 
			 * @param path The path to the folder where the plot will be saved. 
			 */
			PlotDistance(const hist::ScatteringHistogram& data, std::string path);

			/**
			 * @brief Destructor. 
			 */
			~PlotDistance() override;
		
			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
			 */
			static void quick_plot(const hist::ScatteringHistogram& data, std::string path);
	};
}