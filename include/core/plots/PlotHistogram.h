// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/HistFwd.h>
#include <plots/Plot.h>
#include <utility/observer_ptr.h>

namespace ausaxs::plots {
	/**
	 * @brief Plot a specific \class Histogram object.
	 */
	class PlotHistogram : public Plot {
		public:
			PlotHistogram();

			~PlotHistogram() override;

			PlotHistogram(const hist::Histogram& h, const plots::PlotOptions& options);

			PlotHistogram& plot(const hist::Histogram& hist, const plots::PlotOptions& options);

			/**
			 * @brief Plot and save the input dataset and the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a histogram. 
			 */
			static void quick_plot(const hist::Histogram& hist, const plots::PlotOptions& options, const io::File& path);
	};
}