#pragma once

#include <plots/PlotOptions.h>
#include <dataset/Dataset.h>
#include <dataset/Multiset.h>
#include <hist/Histogram.h>

#include <string.h>
#include <memory.h>
#include <sstream>

namespace plots {
	/**
	 * @brief \class Plot.
	 * 
	 * Virtual super-class for all plotter objects. 
	 */
	class Plot {
		public: 
			/**
			 * @brief Default constructor.
			 */
			Plot() = default;

			/**
			 * @brief Destructor.
			 */
			virtual ~Plot() = default;

			/**
			 * @brief Write this plot to a given destination. 
			 * @param folder Path to the folder where this plot will be saved. 
			 */
			void save(std::string folder) const;

			std::stringstream ss;
	};
}