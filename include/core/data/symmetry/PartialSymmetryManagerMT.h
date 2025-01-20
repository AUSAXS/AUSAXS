#pragma once

#include <hist/distribution/GenericDistribution1D.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/detail/CompactCoordinates.h>
#include <data/symmetry/detail/SymmetryHelpers.h>

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
			observer_ptr<const data::Molecule> protein;													// the molecule we are calculating the histogram for
            detail::MasterHistogram<use_weighted_distribution> master;									// the current total histogram
			std::vector<symmetry::detail::CompactCoordinateSymmetries> coords;							// a compact representation of the relevant data from the managed bodies
			container::Container2D<detail::PartialHistogram<use_weighted_distribution>> partials_aa; 	// the partial histograms
			container::Container1D<detail::HydrationHistogram<use_weighted_distribution>> partials_aw;	// the partial hydration-atom histograms
			detail::HydrationHistogram<use_weighted_distribution> partials_ww;               			// the partial histogram for the hydration layer

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
			void calc_aa_self(calculator_t calculator, int index);

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_ww_self(calculator_t calculator);

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aa(calculator_t calculator, int n, int m);

			/**
			 * @brief Calculate the atom-atom distances between the symmetric duplicate @a i and @a j of bodies @a n and @a m. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_ss(calculator_t calculator, int n, int m);

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aw(calculator_t calculator, int index);

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
			void update_compact_representation_water(int index);
	};
}