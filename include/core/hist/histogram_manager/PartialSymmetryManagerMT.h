// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/GenericDistribution1D.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/histogram_manager/detail/SymmetryDetailFwd.h>

#include <memory>
#include <mutex>

namespace ausaxs::hist {
	/**
	 * @brief A multi-threaded smart distance calculator which efficiently calculates the simple distance histogram. 
	 */
    template<bool use_weighted_distribution> 
	class PartialSymmetryManagerMT : public IPartialHistogramManager {
		using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
		using calculator_t = observer_ptr<distance_calculator::SimpleCalculator<use_weighted_distribution>>;

		template<typename T> using BodyIndexer2D = typename container::Container2D<T>;
		template<typename T> using BodyIndexer1D = typename container::Container1D<T>;

		// 2D symmetry indexer to be stored within a BodyIndexer2D
		template<typename T> struct SymmetryIndexer2D {
			SymmetryIndexer2D() = default;
			SymmetryIndexer2D(int size, T&& value) : data(size, std::vector<T>(size, std::forward<T>(value))) {}
			SymmetryIndexer2D(int size_x, int size_y, T&& value) : data(size_x, std::vector<T>(size_y, std::forward<T>(value))) {}
			T& index(int isym1, int isym2) {return data[isym1][isym2];}
			std::vector<std::vector<T>> data;
		}; 

		// 1D symmetry indexer to be stored within a BodyIndexer1D
		template<typename T> struct SymmetryIndexer1D {
			SymmetryIndexer1D() = default;
			SymmetryIndexer1D(int size, T&& value) : data(size, std::forward<T>(value)) {}
			template<typename ...Arg> SymmetryIndexer1D(Arg&&... args) : data(std::forward<Arg>(args)...) {}
			T& index(int isym) {return data[isym];}
			std::vector<T> data;
		};

		public:
			PartialSymmetryManagerMT(observer_ptr<const data::Molecule> protein);
			virtual ~PartialSymmetryManagerMT() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			std::unique_ptr<ICompositeDistanceHistogram> calculate_all() override;

		private:
			observer_ptr<const data::Molecule> protein;					// the molecule we are calculating the histogram for
            detail::MasterHistogram<use_weighted_distribution> master;	// the current total histogram
			std::vector<symmetry::detail::BodySymmetryData> coords;		// a compact representation of the relevant data from the managed bodies
			hist::detail::CompactCoordinates coords_w;					// a compact representation of the relevant data from the hydration layer
			std::unordered_map<int, int> res_self_index_map;			// a map to keep track of result indexes in the self-correlation results
			std::unordered_map<int, int> res_cross_index_map;			// a map to keep track of result indexes in the cross-correlation results

			// partial histograms - the types are quite complex since we must track both bodies and symmetries
			BodyIndexer2D<SymmetryIndexer2D<detail::PartialHistogram<use_weighted_distribution>>> partials_aa;
			BodyIndexer1D<SymmetryIndexer1D<detail::HydrationHistogram<use_weighted_distribution>>> partials_aw;
			detail::HydrationHistogram<use_weighted_distribution> partials_ww;
			std::mutex master_hist_mutex;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			template<bool hydration_enabled>
			std::unique_ptr<DistanceHistogram> _calculate();

			/**
			 * @brief Initialize the storage spaces of this object. 
			 */
			void initialize();

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
			void calc_aa_self(calculator_t calculator, int ibody) const;

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_ww(calculator_t calculator) const;

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aa(calculator_t calculator, int ibody1, int isym1, int ibody2, int isym2) const;

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 *
			 * @param ibody The index of the body to calculate the self-correlation for.
			 * @param isym The index of the symmetry to calculate the self-correlation for. Index 0 is the main body.
			 */
			void calc_aw(calculator_t calculator, int ibody, int isym) const;

			/**
			 * @brief Combine the atom-atom correlation of two bodies into the master histogram.
			 *
			 * @param ibody1 The index of the first body to combine the correlation for.
			 * @param isym1 The index of the symmetry of the first body to combine the correlation for. Index 0 is the main body.
			 * @param ibody2 The index of the second body to combine the correlation for.
			 * @param isym2 The index of the symmetry of the second body to combine the correlation for. Index 0 is the main body.
			 */
			void combine_aa(int ibody1, int isym1, int ibody2, int isym2, GenericDistribution1D_t&& res);

			/**
			 * @brief Combine the atom-hydration correlation of a body symmetry into the master histogram.
			 *
			 * @param ibody The index of the body to combine the correlation for.
			 * @param isym The index of the symmetry to combine the correlation for. Index 0 is the main body.
			 */
			void combine_aw(int ibody, int isym, GenericDistribution1D_t&& res);

			/**
			 * @brief Combine the hydration-hydration correlation into the master histogram.
			 */
			void combine_ww(GenericDistribution1D_t&& res);

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