// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/histogram_manager/IPartialHistogramManager.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/detail/SimpleExvModel.h>
#include <container/Container1D.h>
#include <container/Container2D.h>

namespace ausaxs::hist {
	/**
	 * The basic idea is that we have a bunch of partial histograms (contained in @a partials), which combined represents the total scattering histogram. 
	 * As an example, if we had 4 bodies, it would look something like this:
	 * 4       x
	 * 3     x
	 * 2   x
	 * 1 x
	 *   1 2 3 4
	 * The self-correlation partials are marked with an 'x'. They are constant and are thus precalculated when this class is initialized. 
	 * The uaaer and lower triangle are symmetric, and we can thus just calculate one of them and double the result. After all partials are initially
	 * generated, this class recalculates them whenever a body has changed. If body 2 is moved, the partials (1, 2), (2, 3), and (2, 4) must be recalculated. 
	 * 
	 * This is further complicated by the presence of the hydration layer. Since this does not belong to any individual body, it can be viewed as 
	 * a simple extension to the above example, so we now have {1, 2, 3, 4, H}. 
	 */

	/**
	 * @brief A single-threaded smart distance calculator which efficiently calculates the simple distance histogram.
	 */
    template<bool weighted_bins, bool variable_bin_width> 
	class PartialHistogramManager : public IPartialHistogramManager {
	    using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
		public:
			PartialHistogramManager(observer_ptr<const data::Molecule> protein); 
			virtual ~PartialHistogramManager() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			virtual std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			virtual std::unique_ptr<ICompositeDistanceHistogram> calculate_all() override;

		protected:
			observer_ptr<const data::Molecule> protein;										// the molecule we are calculating the histogram for
            detail::MasterHistogram<weighted_bins> master;									// the current total histogram
            std::vector<detail::CompactCoordinates<variable_bin_width>> coords_a;			// a compact representation of the relevant data from the managed bodies
            detail::CompactCoordinates<variable_bin_width> coords_w;                		// a compact representation of the hydration data
			container::Container2D<detail::PartialHistogram<weighted_bins>> partials_aa; 	// the partial histograms
			container::Container1D<detail::HydrationHistogram<weighted_bins>> partials_aw;	// the partial hydration-atom histograms
			detail::HydrationHistogram<weighted_bins> partials_ww;               			// the partial histogram for the hydration layer

		private:
			/**
			 * @brief Initialize this object. The internal distances between atoms in each body is constant and cannot change. 
			 *        They are unaffected by both rotations and translations, and so we precalculate them. 
			 */
			void initialize();

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 */
			void calc_aa(unsigned int n, unsigned int m);

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 */
			void calc_aw(unsigned int index);

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 */
			void calc_ww();

			/**
			 * @brief Calculate the self-correlation of a body. 
			 */
			void calc_self_correlation(unsigned int index);
	};
}