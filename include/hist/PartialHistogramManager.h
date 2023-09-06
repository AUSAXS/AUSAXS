#pragma once

#include <hist/HistogramManager.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/detail/CompactCoordinates.h>
#include <utility/Container1D.h>
#include <utility/Container2D.h>

namespace hist {
	/**
	 * The basic idea is that we have a bunch of partial histograms (contained in @a partials), which combined represents the total scattering histogram. 
	 * As an example, if we had 4 bodies, it would look something like this:
	 * 4       x
	 * 3     x
	 * 2   x
	 * 1 x
	 *   1 2 3 4
	 * The self-correlation partials are marked with an 'x'. They are constant and are thus precalculated when this class is initialized. 
	 * The upper and lower triangle are symmetric, and we can thus just calculate one of them and double the result. After all partials are initially
	 * generated, this class recalculates them whenever a body has changed. If body 2 is moved, the partials (1, 2), (2, 3), and (2, 4) must be recalculated. 
	 * 
	 * This is further complicated by the presence of the hydration layer. Since this does not belong to any individual body, it can be viewed as 
	 * a simple extension to the above example, so we now have {1, 2, 3, 4, H}. 
	 */

	/**
	 * @brief A single-threaded smart distance calculator which efficiently calculates the distance histogram.
	 */
	class PartialHistogramManager : public HistogramManager {
		public:
			PartialHistogramManager(Protein* protein); 

			virtual ~PartialHistogramManager() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			virtual Histogram calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			virtual ScatteringHistogram calculate_all() override;

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
			std::vector<detail::CompactCoordinates> coords_p;   // a compact representation of the relevant data from the managed bodies
			detail::CompactCoordinates coords_h;                // a compact representation of the hydration data

			// histogram data
			detail::MasterHistogram master;                       	// the current total histogram
			Container2D<detail::PartialHistogram> partials_pp; 		// the partial histograms
			Container1D<detail::HydrationHistogram> partials_hp;	// the partial hydration-atom histograms
			detail::HydrationHistogram partials_hh;               	// the partial histogram for the hydration layer

			/**
			 * @brief Initialize this object. The internal distances between atoms in each body is constant and cannot change. 
			 *        They are unaffected by both rotations and translations, and so we precalculate them. 
			 */
			void initialize();

			/**
			 * @brief Calculate the atom-atom distances between body @a index and all others. 
			 */
			void calc_pp(unsigned int index);

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 */
			void calc_pp(unsigned int n, unsigned int m);

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 */
			void calc_hp(unsigned int index);

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 */
			void calc_hh();

			/**
			 * @brief Calculate the self-correlation of a body. 
			 */
			void calc_self_correlation(unsigned int index);
	};
}