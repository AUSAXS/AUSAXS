#pragma once

#include <plots/Plot.h>
#include <hist/Histogram.h>

namespace plots {
	class PlotHistogram : public Plot {
		public:
			/**
			 * @brief Copy constructor.
			 * 
			 * @param h The Histogram to be plotted. 
			 */
			PlotHistogram(const hist::Histogram& h);

			/**
			 * @brief Destructor.
			 */
			~PlotHistogram() override;

			void plot(const hist::Histogram& hist);

			/**
			 * @brief Plot and save the input dataset and the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a histogram. 
			 */
			static void quick_plot(const hist::Histogram& hist, const io::File& path);
	};
}