#pragma once

#include <hist/detail/BodyTracker.h>
#include <data/DataFwd.h>

#include <vector>
#include <memory>

namespace hist {
	class DistanceHistogram;
	class CompositeDistanceHistogram;

	/**
	 * @brief A histogram manager which calculates the distance histogram in a slow but simple way. 
	 * 		  This class is only intended for testing and inheritance. Use the PartialHistogramManagerMT class for production. 
	 */
	class HistogramManager : public hist::BodyTracker {
		public:
			HistogramManager(const data::Molecule* const protein); 

			HistogramManager(const HistogramManager& hm); 

			virtual ~HistogramManager();

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			virtual std::unique_ptr<DistanceHistogram> calculate();

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			virtual std::unique_ptr<CompositeDistanceHistogram> calculate_all();

		protected:
			const data::Molecule* const protein;	// pointer to the parent Protein
    };
}