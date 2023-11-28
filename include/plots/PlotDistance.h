#pragma once

#include <plots/Plot.h>
#include <hist/HistFwd.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

namespace plots {
	/**
	 * @brief Plot the distance histogram for a protein. 
	 */
	class PlotDistance : public Plot {
		public:
			/**
			 * @brief Constructor.
			 * 
			 * @param data The ScatteringHistogram which will be plotted. 
			 * @param path The path to the folder where the plot will be saved. 
			 */
			PlotDistance(observer_ptr<hist::ICompositeDistanceHistogram> data, const io::File& path);

			/**
			 * @brief Constructor.
			 * 
			 * @param data The protein whose distances will be plotted. 
			 * @param path The path to the folder where the plot will be saved. 
			 */
			PlotDistance(observer_ptr<data::Molecule> protein, const io::File& path);

			/**
			 * @brief Destructor. 
			 */
			~PlotDistance() override;
		
			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
			 */
			static void quick_plot(observer_ptr<hist::ICompositeDistanceHistogram> data, const io::File& path);
		
			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
			 */
			static void quick_plot(observer_ptr<data::Molecule> protein, const io::File& path);
	};
}