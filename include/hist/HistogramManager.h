#pragma once

// forwards declaration
class Protein;

#include <data/StateManager.h>
#include <hist/ScatteringHistogram.h>
#include <hist/Histogram.h>
#include <hist/detail/BodyTracker.h>

#include <vector>

namespace hist {
	/**
	 * @brief A histogram manager which calculates the distance histogram in a slow but simple way. 
	 * 		  This class is only intended for testing and inheritance. Use the PartialHistogramManagerMT class for production. 
	 */
	class HistogramManager : public hist::BodyTracker {
		public:
			HistogramManager(Protein* protein); 

			virtual ~HistogramManager();

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			virtual Histogram calculate();

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			virtual ScatteringHistogram calculate_all();

			enum class Type {
				Simple,
				SimpleMT,
				Partial,
				PartialMT
			};

		protected:
			Protein* protein;	// pointer to the parent Protein
    };
}