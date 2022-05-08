#pragma once

#include <plots/Plot.h>
#include <histogram/Histogram.h>

#include <memory>
#include <string>

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

            /**
             * @brief Save this image at the given location in the specified format. 
             * 
             * @param path The path & format of the image. 
             */
			void save(std::string path) const override;

			/**
			 * @brief Plot and save the input dataset and the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a histogram. 
			 */
			static void quick_plot(const hist::Histogram& hist, std::string path);

		private:
			std::shared_ptr<TCanvas> canvas;

			void initial_plot(const hist::Histogram& hist);

			void prepare_canvas();
		};
}