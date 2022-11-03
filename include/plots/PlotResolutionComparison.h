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
			PlotResolutionComparison(Multiset d);

			/**
			 * @brief Destructor.
			 */
			~PlotResolutionComparison() override;

			/**
			 * @brief Plot and save the input dataset and the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single dataset. 
			 */
			static void quick_plot(const Multiset& data, std::string path);
	};
}