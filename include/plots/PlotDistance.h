#pragma once

#include <plots/Plot.h>
#include <histogram/ScatteringHistogram.h>

namespace plots {
	/**
	 * @brief \class PlotDistance.
	 * 
	 * Plots a histogram of all distances. 
	 */
	class PlotDistance : public Plot {
		public:
		/**
		 * @brief Constructor.
		 * @param d The ScatteringHistogram which will be plotted. 
		 */
		PlotDistance(const histogram::ScatteringHistogram& d) : d(d) {}

		/**
		 * @brief Destructor. 
		 */
		~PlotDistance() override;

		/**
		 * @brief Save this plot at the given location.
		 * @param path Path to where this plot will be saved. Note that it will attempt to save in the specified format. 
		 */
		void save(std::string path) const override;

		private: 
		const histogram::ScatteringHistogram d; // The ScatteringHistogram backing this object. 
	};
}