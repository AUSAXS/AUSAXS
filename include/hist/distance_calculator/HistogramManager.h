#pragma once

#include <hist/distance_calculator/IHistogramManager.h>
#include <data/DataFwd.h>
#include <utility/view_ptr.h>

#include <vector>
#include <memory>

namespace hist {
	/**
	 * @brief A histogram manager which calculates the distance histogram in a slow but simple way. 
	 * 		  This class is only intended for testing and inheritance. Use the PartialHistogramManagerMT class for production. 
	 */
	template<bool use_weighted_distribution>
	class HistogramManager : public IHistogramManager {
		public:
			HistogramManager(view_ptr<const data::Molecule> protein); 

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
			view_ptr<const data::Molecule> protein; // pointer to the parent Protein
    };
}