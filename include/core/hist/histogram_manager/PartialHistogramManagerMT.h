// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/GenericDistribution1D.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/detail/CompactCoordinates.h>

#include <memory>
#include <mutex>

namespace ausaxs::hist {
	/**
	 * @brief A multi-threaded smart distance calculator which efficiently calculates the simple distance histogram. 
	 */
    template<bool weighted_bins, bool variable_bin_width> 
	class PartialHistogramManagerMT : public PartialHistogramManager<weighted_bins, variable_bin_width> {
		public:
			PartialHistogramManagerMT(observer_ptr<const data::Molecule> protein);
			virtual ~PartialHistogramManagerMT() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			std::unique_ptr<ICompositeDistanceHistogram> calculate_all() override;

		private:
		    using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
			using calculator_t = observer_ptr<distance_calculator::SimpleCalculator<weighted_bins, variable_bin_width>>;
			struct { // cache for early return
				GenericDistribution1D_t p_aa;
				GenericDistribution1D_t p_aw;
				GenericDistribution1D_t p_ww;
				GenericDistribution1D_t p_tot;
			} cache;
			std::mutex master_hist_mutex;

			/**
			 * @brief Initialize this object. The internal distances between atoms in each body is constant and cannot change. 
			 *        They are unaffected by both rotations and translations, and so we precalculate them. 
			 */
			void initialize(calculator_t calculator);

			/**
			 * @brief Calculate the self-correlation of a body.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_self_correlation(calculator_t calculator, int index);

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aa(calculator_t calculator, int n, int m);

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aw(calculator_t calculator, int index);

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_ww(calculator_t calculator);

			void combine_self_correlation(int index, GenericDistribution1D_t&&);

			void combine_aa(int n, int m, GenericDistribution1D_t&&);

			void combine_aw(int index, GenericDistribution1D_t&&);

			void combine_ww(GenericDistribution1D_t&&);

			/**
			 * @brief Update the compact representation of the coordinates of body @a index.
			 * 
			 * @param index The index of the body to update.
			 */
			void update_compact_representation_body(int index);

			/**
			 * @brief Update the compact representation of the coordinates of the hydration layer.
			 */
			void update_compact_representation_water();
	};
}