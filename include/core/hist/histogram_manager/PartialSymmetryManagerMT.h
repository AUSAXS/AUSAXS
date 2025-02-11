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
			 * @brief Initialize this object. The internal distances between atoms in each body is constant and cannot change. 
			 *        They are unaffected by both rotations and translations, and so we precalculate them. 
			 */
			void initialize(calculator_t calculator);

			/**
			 * @brief Calculate the self-correlation of a body.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aa_self(calculator_t calculator, int ibody);

			void calc_aa_self(calculator_t calculator, int ibody, int isym);

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
			 */
			void calc_aw(calculator_t calculator, int ibody, int isym1);

			void combine_aa_self(int index, GenericDistribution1D_t&&);

			void combine_aa(int ibody1, int isym1, int ibody2, int isym2, GenericDistribution1D_t&&);

			void combine_aw(int ibody, int isym, GenericDistribution1D_t&&);

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