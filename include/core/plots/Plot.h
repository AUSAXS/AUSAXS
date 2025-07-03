// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>
#include <plots/PlotOptions.h>

#include <sstream>

namespace ausaxs::plots {
	/**
	 * @brief Virtual super-class for all plotter objects. 
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
			 * 
			 * @param path Path to the folder where this plot will be saved. 
			 */
			void save(const io::File& path) const;

			std::stringstream ss;
	};
}