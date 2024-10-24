#pragma once

#include <plots/Plot.h>
#include <utility/observer_ptr.h>
#include <hist/HistFwd.h>

namespace ausaxs::plots {
	/**
	 * @brief Plot a specific \class Histogram object.
	 */
	class PlotHistogram : public Plot {
		public:
			PlotHistogram();

			virtual ~PlotHistogram();

			PlotHistogram(const hist::Histogram& h, const plots::PlotOptions& options);

			PlotHistogram& plot(const hist::Histogram& hist, const plots::PlotOptions& options);

			/**
			 * @brief Plot and save the input dataset and the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a histogram. 
			 */
			static void quick_plot(const hist::Histogram& hist, const plots::PlotOptions& options, const io::File& path);
	};
}