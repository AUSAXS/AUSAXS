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
	template<bool use_weighted_distribution>
	class HistogramManager : public IHistogramManager {
		public:
			HistogramManager(observer_ptr<const data::Molecule> protein); 
			virtual ~HistogramManager();

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