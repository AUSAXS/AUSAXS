// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/distance_calculator/DistanceCalculatorFwd.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/histogram_manager/PartialHistogramManager.h>

#include <memory>
#include <vector>

namespace ausaxs::hist {
	/**
	 * @brief A multi-threaded smart distance calculator which efficiently calculates the simple distance histogram. 
	 */
    template<bool weighted_bins, bool variable_bin_width> 
	// NOLINTNEXTLINE - the destructor is virtual through the dependent base, which the check cannot see on the template pattern
	class PartialHistogramManagerMT : public PartialHistogramManager<weighted_bins, variable_bin_width> {
		public:
			PartialHistogramManagerMT(observer_ptr<const data::Molecule> protein);
			~PartialHistogramManagerMT() override;

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
			using calculator_t = observer_ptr<distance_calculator::Calculator<weighted_bins, variable_bin_width>>;
			struct { // cache for early return
				GenericDistribution1D_t p_aa;
				GenericDistribution1D_t p_aw;
				GenericDistribution1D_t p_ww;
				GenericDistribution1D_t p_tot;
			} cache;
			std::unique_ptr<distance_calculator::HistogramStore<weighted_bins>> store;
			std::vector<int> recalculated; // the rows queued for recalculation in the current run, see recalculate()

			int handle_aa(int n, int m) const; // m <= n
			int handle_aw(int index) const;
			int handle_ww() const;

			/**
			 * @brief Initialize the master histogram and the storage of the partial histograms.
			 */
			void initialize(int bin_count);

			/**
			 * @brief Take the partial histogram @a h out of the master histogram before it is recalculated.
			 *        Must be called before the first calculation into @a h is queued; the new contents are added back once
			 *        the calculator has run.
			 */
			void recalculate(int h);

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