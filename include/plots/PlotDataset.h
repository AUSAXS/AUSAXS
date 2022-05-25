#pragma once

#include <plots/Plot.h>
#include <utility/Dataset.h>
#include <utility/Multiset.h>

#include <memory>
#include <string>

namespace plots {
	class PlotDataset : public Plot {
		public:
			/**
			 * @brief Constructor.
			 * 
			 * @param d The Dataset to be plotted. 
			 */
			PlotDataset(const Dataset& d);

			/**
			 * @brief Constructor.
			 * 
			 * @param d The Multiset to be plotted. 
			 */
			PlotDataset(const Multiset& d);

			/**
			 * @brief Destructor.
			 */
			~PlotDataset() override;

			/**
			 * @brief Plot an additional Dataset. 
			 */
			void plot(const Dataset& data);

            /**
             * @brief Save this image at the given location in the specified format. 
             * 
             * @param path The path & format of the image. 
             */
			void save(std::string path) const override;

			/**
			 * @brief Plot and save the input dataset at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
			 */
			static void quick_plot(const Dataset& data, std::string path);

		private:
			std::shared_ptr<TCanvas> canvas;

			void initial_plot(const Dataset& data);

			void prepare_canvas();
		};
}