#pragma once

#include <plots/Plot.h>
#include <dataset/DatasetFwd.h>

namespace ausaxs::plots {
	class PlotOptions;

	/**
	 * @brief Plot the contents of a specific \class Dataset object.
	 * 		  Multiple datasets can be plotted by chaining calls to the plot() method.
	 */
	class PlotDataset : public Plot {
		public:
			/**
			 * @brief Default constructor.
			 */
			PlotDataset() = default;

			/**
			 * @brief Constructor.
			 */
			PlotDataset(const Dataset& d, const plots::PlotOptions& options);

			/**
			 * @brief Destructor.
			 */
			~PlotDataset() override;

			/**
			 * @brief Plot an additional Dataset. 
			 */
			PlotDataset& plot(const Dataset& data, const plots::PlotOptions& options);

			/**
			 * @brief Plot a vertical line at the specified x coordinate.
			 */
			PlotDataset& vline(double x, const plots::PlotOptions& options);

			/**
			 * @brief Plot a horizontal line at the specified y coordinate.
			 */
			PlotDataset& hline(double y, const plots::PlotOptions& options);

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
			 */
			static void quick_plot(const Dataset& data, const plots::PlotOptions& options, const io::File& path);
		};
}