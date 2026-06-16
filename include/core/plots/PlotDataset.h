// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

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
			PlotDataset() = default;
			PlotDataset(const Dataset& d, const plots::PlotOptions& options);
			~PlotDataset() override;

			/**
			 * @brief Plot an additional Dataset. 
			 */
			PlotDataset& plot(const Dataset& data, const plots::PlotOptions& options);

			/**
			 * @brief Plot dataset of the form q | I | Ierr | model as a two-panel data + residuals plot
			 *        This does not support chaining since .plot files currently cannot specify which panel to draw on.
			 */
			void plot_residuals(const Dataset& data, const plots::PlotOptions& options);

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