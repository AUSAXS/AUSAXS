#pragma once

// forwards declaration
class Protein;

#include <data/Atom.h>
#include <data/Body.h>
#include <data/StateManager.h>
#include <hist/ScatteringHistogram.h>
#include <hist/Histogram.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/detail/CompactCoordinates.h>

#include <vector>

namespace hist {
	/**
	 * @brief A histogram manager which calculates the distance histogram in a slow but simple way. 
	 * 		  This class is only intended for testing and inheritance. Use the PartialHistogramManagerMT class for production. 
	 */
	class HistogramManager {
		public:
			HistogramManager(Protein* protein); 

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			virtual Histogram calculate();

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			virtual ScatteringHistogram calculate_all();

			/**
			 * @brief Get a signalling object for signalling a change of state. 
			 *        Each body is supposed to hold one of these, and trigger it when they change state. 
			 */
			std::shared_ptr<StateManager::BoundSignaller> get_probe(unsigned int i);

			/**
			 * @brief Signal that the hydration layer was modified. 
			 *        This is supposed to be used only by the Protein class, which has direct access to this object. Thus a signalling object is unnecessary. 
			 */
			void signal_modified_hydration_layer();

			const StateManager& get_state_manager() const;

			StateManager& get_state_manager();

		protected:
			const unsigned int size;                            // number of managed bodies
			StateManager statemanager;                    		// a helper which keeps track of state changes in each body
			std::vector<detail::CompactCoordinates> coords_p;   // a compact representation of the relevant data from the managed bodies
			detail::CompactCoordinates coords_h;                // a compact representation of the hydration data
			Protein* protein;                             		// pointer to the parent Protein

			// histogram data
			detail::MasterHistogram master;                       			// the current total histogram
			std::vector<std::vector<detail::PartialHistogram>> partials_pp; // the partial histograms
			std::vector<detail::HydrationHistogram> partials_hp;       		// the partial hydration-atom histograms
			detail::HydrationHistogram partials_hh;               			// the partial histogram for the hydration layer
    };
}