// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/distance_calculator/DistanceCalculatorFwd.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/detail/SymmetryDetailFwd.h>

#include <cassert>
#include <memory>
#include <vector>

namespace ausaxs::hist {
	/**
	 * @brief A multi-threaded smart distance calculator which efficiently calculates the simple distance histogram. 
	 */
    template<bool weighted_bins, bool variable_bin_width> 
	class PartialSymmetryManagerMT : public IPartialHistogramManager {
		public:
			PartialSymmetryManagerMT(observer_ptr<const data::Molecule> protein);
			~PartialSymmetryManagerMT() override;

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

			observer_ptr<const data::Molecule> protein;									// the molecule we are calculating the histogram for
            detail::MasterHistogram<weighted_bins> master;								// the current total histogram
			std::vector<symmetry::detail::BodySymmetryData<variable_bin_width>> coords;	// a compact representation of the relevant data from the managed bodies
			hist::detail::CompactCoordinates<variable_bin_width> coords_w;				// a compact representation of the relevant data from the hydration layer
			std::unique_ptr<distance_calculator::HistogramStore<weighted_bins>> store;
			std::vector<int> recalculated; // the results queued for recalculation in the current run, see recalculate()

			/**
			 * @brief The store ids of the atom-atom results, one for each pair of symmetries of each body pair.
			 *        The symmetry index isym runs over the main body (0) and its symmetries (1..).
			 *        Only the lower triangle is calculated: ibody2 <= ibody1, and within a body isym2 < isym1 or (0, 0).
			 *        The other pairs get no result, so they cost no histogram.
			 */
			class SymmetryPairIds {
				public:
					SymmetryPairIds() = default;

					/**
					 * @brief Lay out the pairs of bodies with @a sym_counts symmetry indices each, and give each calculated pair the id returned by @a allocate.
					 */
					template<typename Allocate>
					SymmetryPairIds(std::vector<int> sym_counts, Allocate&& allocate) : sym_counts(std::move(sym_counts)) {
						int n = static_cast<int>(this->sym_counts.size());
						for (int ibody1 = 0; ibody1 < n; ++ibody1) {
							for (int ibody2 = 0; ibody2 <= ibody1; ++ibody2) {
								offsets.push_back(static_cast<int>(ids.size()));
								for (int isym1 = 0; isym1 < this->sym_counts[ibody1]; ++isym1) {
									for (int isym2 = 0; isym2 < this->sym_counts[ibody2]; ++isym2) {
										ids.push_back(calculated(ibody1, isym1, ibody2, isym2) ? allocate() : -1);
									}
								}
							}
						}
					}

					int id(int ibody1, int isym1, int ibody2, int isym2) const {
						assert(0 <= ibody2 && ibody2 <= ibody1 && ibody1 < static_cast<int>(sym_counts.size()) && "SymmetryPairIds::id: expected a body pair in the lower triangle");
						assert(0 <= isym1 && isym1 < sym_counts[ibody1] && 0 <= isym2 && isym2 < sym_counts[ibody2] && "SymmetryPairIds::id: symmetry index out of range; symmetries may not be added after the first calculation");
						int id = ids[offsets[ibody1*(ibody1+1)/2 + ibody2] + isym1*sym_counts[ibody2] + isym2];
						assert(id != -1 && "SymmetryPairIds::id: expected a symmetry pair in the lower triangle");
						return id;
					}

					/**
					 * @brief Call @a f with the id of every calculated pair.
					 */
					template<typename F>
					void for_each_id(F&& f) const {
						for (int id : ids) {if (id != -1) {f(id);}}
					}

				private:
					static bool calculated(int ibody1, int isym1, int ibody2, int isym2) {
						return ibody1 != ibody2 || isym2 < isym1 || (isym1 == 0 && isym2 == 0);
					}

					std::vector<int> sym_counts; // [ibody]
					std::vector<int> offsets;    // [ibody1*(ibody1+1)/2 + ibody2], the start of the pair's block in ids
					std::vector<int> ids;        // each block is [isym1][isym2]; -1 for a pair that is not calculated
			};

			// the result ids in the store, fixed by initialize()
			SymmetryPairIds aa;
			std::vector<std::vector<int>> aw; // [ibody][isym]
			int ww = -1;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			template<bool hydration_enabled>
			std::unique_ptr<DistanceHistogram> _calculate();

			/**
			 * @brief Determine the number of bins, discarding everything calculated so far if the structure outgrew them.
			 */
			int prepare_axis();

			/**
			 * @brief Initialize the master histogram, the layout of the partial histograms, and their storage.
			 *        The number of symmetries of each body is fixed from here on.
			 */
			void initialize(int bin_count);

			/**
			 * @brief Take the partial histogram @a id out of the master histogram before it is recalculated.
			 */
			void recalculate(int id);

			/**
			 * @brief Expand the modification flags for shared reference symmetries.
			 */
			void propagate_reference_symmetry_modifications(
				const std::vector<bool>& externally_modified,
				const std::vector<bool>& internally_modified,
				std::vector<std::vector<bool>>& symmetry_modified
			) const;

			/**
			 * @brief Calculate the self-correlation of a body. 
			 *		  This includes: 
			 *		      1. The self-correlation of the main body.
			 *		      2. The self-correlation of each symmetry of the body.
			 *		  No internal cross terms are calculated here.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 *
			 * @param ibody The index of the body to calculate the self-correlation for.
			 */
			void calc_aa_self(calculator_t calculator, int ibody);

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_ww(calculator_t calculator);

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aa(calculator_t calculator, int ibody1, int isym1, int ibody2, int isym2);

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 *
			 * @param ibody The index of the body to calculate the self-correlation for.
			 * @param isym The index of the symmetry to calculate the self-correlation for. Index 0 is the main body.
			 */
			void calc_aw(calculator_t calculator, int ibody, int isym);

			/**
			 * @brief Update the compact representation of the coordinates of body @a index.
			 * 
			 * @param index The index of the body to update.
			 */
			void update_compact_representation_body(int ibody);

			/**
			 * @brief Update the compact representation of the coordinates of body @a index.
			 * 
			 * @param ibody The index of the body to update.
			 * @param isym The index of the symmetry to update. Index 0 is the main body.
			 */
			void update_compact_representation_symmetry(int ibody, int isym);

			/**
			 * @brief Update the compact representation of the coordinates of the hydration layer.
			 */
			void update_compact_representation_water();
	};
}