// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/distance_calculator/DistanceCalculatorFwd.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/histogram_manager/IPartialHistogramManager.h>
#include <settings/ExvSettings.h>

#include <type_traits>

#include <memory>
#include <vector>

namespace ausaxs::hist {
	/**
	 * @brief Common machinery of the multithreaded partial histogram managers, which only recalculate the parts of the histogram
	 *        changed between each call.
	 *
	 * This is independent of the single-threaded PartialHistogramManager, which is kept simple as a reference implementation.
	 *
	 * @tparam form_factors Whether the atoms are resolved by form factor, see HistogramManagerMTBase.
	 */
    template<bool weighted_bins, bool form_factors> 
	class PartialHistogramManagerMTBase : public IPartialHistogramManager {
		public:
			/**
			 * @param exv_method The excluded volume model the form factor-resolved result is built for; see detail::make_histogram.
			 */
			PartialHistogramManagerMTBase(observer_ptr<const data::Molecule> protein, settings::exv::ExvMethod exv_method);
			~PartialHistogramManagerMTBase() override;

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
			using calculator_t = observer_ptr<distance_calculator::Calculator<weighted_bins, form_factors>>;
			using AtomicCoordinates = std::conditional_t<form_factors, std::vector<hist::detail::CompactCoordinates>, hist::detail::CompactCoordinates>;

			observer_ptr<const data::Molecule> protein;		// the molecule we are calculating the histogram for
			settings::exv::ExvMethod exv_method;			// the excluded volume model of the form factor-resolved result
			detail::MasterHistogram<weighted_bins> master;	// the current total histogram
			std::vector<AtomicCoordinates> coords_a;		// a compact representation of the atoms of each body; with form factors split by type
			hist::detail::CompactCoordinates coords_w;		// a compact representation of the hydration layer
			GenericDistribution1D_t cached_p_tot; // the total histogram of the last calculation, returned as is while nothing is modified
			std::unique_ptr<distance_calculator::HistogramStore<weighted_bins>> store;
			std::vector<std::vector<int>> aa; // the result ids in the store per body pair [n][m], only calculated for m <= n
			std::vector<int> aw;              // the result ids in the store per body
			int ww = -1;                      // the result id in the store of the hydration layer
			std::vector<int> recalculated;    // the results queued for recalculation in the current run, see recalculate()

			/**
			 * @brief Determine the number of bins, discarding everything calculated so far if the structure outgrew them.
			 */
			int prepare_axis();

			/**
			 * @brief Initialize the master histogram and the storage of the partial histograms.
			 */
			void initialize(int bin_count);

			/**
			 * @brief Take the partial histogram @a id out of the master histogram before it is recalculated.
			 *        The new contents are added back once the calculator has run. With form factors, that is every histogram of its classes.
			 */
			void recalculate(int id);

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

	/**
	 * @brief The partial histogram manager for the simple excluded volume model, where every atom carries its own weight.
	 */
	template<bool weighted_bins>
	// NOLINTNEXTLINE - the destructor is virtual through the dependent base, which the check cannot see on the template pattern
	class PartialHistogramManagerMT : public PartialHistogramManagerMTBase<weighted_bins, false> {
		public:
			explicit PartialHistogramManagerMT(observer_ptr<const data::Molecule> protein) 
				: PartialHistogramManagerMTBase<weighted_bins, false>(protein, settings::exv::ExvMethod::Simple) {}
	};

	/**
	 * @brief The partial histogram manager for the form factor-resolved excluded volume models.
	 */
	template<bool weighted_bins>
	// NOLINTNEXTLINE - the destructor is virtual through the dependent base, which the check cannot see on the template pattern
	class PartialHistogramManagerMTFF : public PartialHistogramManagerMTBase<weighted_bins, true> {
		public:
			explicit PartialHistogramManagerMTFF(observer_ptr<const data::Molecule> protein, settings::exv::ExvMethod exv_method = settings::exv::exv_method) 
				: PartialHistogramManagerMTBase<weighted_bins, true>(protein, exv_method) {}
	};
}