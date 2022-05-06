#pragma once

#include <plots/Plot.h>
#include <histogram/ScatteringHistogram.h>

#include <memory>
#include <string>

namespace plots {
	class PlotDataset : public Plot {
		public:
			/**
			 * @brief Copy constructor.
			 * 
			 * @param d The ScatteringHistogram to be plotted. 
			 */
			PlotDataset(const Dataset& d);

			/**
			 * @brief Destructor.
			 */
			~PlotDataset() override;

			/**
			 * @brief Plot an additional data set as points. 
			 */
			void plot(const Dataset& data);

            /**
             * @brief Save this image at the given location in the specified format. 
             * 
             * @param path The path & format of the image. 
             */
			void save(std::string path) const override;

			/**
			 * @brief Plot and save the input dataset and the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
			 */
			static void quick_plot(const Dataset& data, std::string path);

		private:
			std::shared_ptr<TCanvas> canvas;

			void initial_plot(const Dataset& data);

			void prepare_canvas();
		};
}