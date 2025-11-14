// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/histogram_manager/IHistogramManager.h>
#include <hist/detail/HistDetailFwd.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::hist {
	/**
	 * @brief A single-threaded simple distance calculator. 
	 *
	 * This class does not account for the excluded volume in any way. 
	 * To implicitly include it, subtract the average excluded volume charge from each atom. 
	 *
	 * This class is not meant for production use. 
	 */
	template<bool weighted_bins, bool variable_bin_width>
	class HistogramManager : public IHistogramManager {
		public:
			virtual ~HistogramManager();
			HistogramManager(observer_ptr<const data::Molecule> protein); 

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			virtual std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			virtual std::unique_ptr<ICompositeDistanceHistogram> calculate_all() override;

		protected:
			observer_ptr<const data::Molecule> protein;
    };
}