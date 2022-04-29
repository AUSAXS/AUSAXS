#pragma once

#include <plots/Plot.h>
#include <ScatteringHistogram.h>

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

			void save(std::string path) const override;

		private:
			std::unique_ptr<TCanvas> canvas;

			void initial_plot(const Dataset& data);

			void prepare_canvas();
		};
}