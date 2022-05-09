#pragma once

#include <plots/Plot.h>
#include <plots/PlotDataset.h>
#include <utility/Multiset.h>

#include <memory>
#include <string>

namespace plots {
	class PlotResolutionComparison : public Plot {
		public:
			/**
			 * @brief Constructor.
			 * 
			 * @param d The Multiset to be plotted. 
			 */
			PlotResolutionComparison(Multiset d, int color = kSolar);

			/**
			 * @brief Destructor.
			 */
			~PlotResolutionComparison() override;

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
			static void quick_plot(const Multiset& data, std::string path);

		private:
			std::shared_ptr<TCanvas> canvas;
			std::unique_ptr<PlotDataset> raw;

            /**
             * @brief Plot the first Dataset. 
             */
			void initial_plot(Dataset& data);

			/**
			 * @brief Plot an additional Dataset. 
			 */
			void plot(Dataset& data);

            /**
             * @brief Prepare the canvas for plotting.
             */
			void prepare_canvas();
		};
}