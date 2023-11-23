#pragma once

#include <plots/PlotHistogram.h>
#include <hist/HistFwd.h>
#include <utility/observer_ptr.h>

namespace plots {
	/**
	 * @brief Plot the distance histogram for a protein. 
	 */
	class PlotProfiles : public PlotHistogram {
		public:
			/**
			 * @brief Constructor.
			 * 
			 * @param data The ScatteringHistogram which will be plotted. 
			 * @param path The path to the folder where the plot will be saved. 
			 */
			PlotProfiles(std::observer_ptr<hist::ICompositeDistanceHistogram> data, const io::File& path);

			/**
			 * @brief Destructor. 
			 */
			~PlotProfiles() override;
		
			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
			 */
			static void quick_plot(std::observer_ptr<hist::ICompositeDistanceHistogram> data, const io::File& path);
	};
}