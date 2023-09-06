#pragma once

#include <hist/HistogramManager.h>
#include <hist/detail/MasterHistogram.h>

namespace hist {
	/**
	 * @brief Base class for the partial histogram managers.
	 */
	class PartialHistogramManagerBase : public HistogramManager {
		public:
			PartialHistogramManagerBase(Protein* protein); 

			virtual ~PartialHistogramManagerBase() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			virtual Histogram calculate() override = 0;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			virtual ScatteringHistogram calculate_all() override = 0;

			/**
			 * @brief Get a signalling object for signalling a change of state. 
			 *        Each body is supposed to hold one of these, and trigger it when they change state. 
			 */
			std::shared_ptr<signaller::Signaller> get_probe(unsigned int i);

			/**
			 * @brief Signal that the hydration layer was modified. 
			 *        This is supposed to be used only by the Protein class, which has direct access to this object. Thus a signalling object is unnecessary. 
			 */
			void signal_modified_hydration_layer();

			const StateManager& get_state_manager() const;

			StateManager& get_state_manager();

		protected:
			detail::MasterHistogram master; // the current total histogram
	};
}