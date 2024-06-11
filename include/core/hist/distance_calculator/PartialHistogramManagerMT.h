#pragma once

#include "hist/distribution/GenericDistribution1D.h"
#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/detail/CompactCoordinates.h>
#include <container/ThreadLocalWrapper.h>
#include <container/Container1D.h>
#include <container/Container2D.h>

#include <memory>
#include <mutex>

namespace hist {
	/**
	 * @brief A multi-threaded smart distance calculator which efficiently calculates the simple distance histogram. 
	 */
    template<bool use_weighted_distribution> 
	class PartialHistogramManagerMT : public PartialHistogramManager<use_weighted_distribution> {
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
		    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
			container::ThreadLocalWrapper<container::Container2D<GenericDistribution1D_t>> partials_aa_all;
			container::ThreadLocalWrapper<container::Container1D<GenericDistribution1D_t>> partials_aw_all;
			container::ThreadLocalWrapper<						 GenericDistribution1D_t>  partials_ww_all;
			std::mutex master_hist_mutex;

			/**
			 * @brief Initialize this object. The internal distances between atoms in each body is constant and cannot change. 
			 *        They are unaffected by both rotations and translations, and so we precalculate them. 
			 */
			void initialize();

			/**
			 * @brief Calculate the atom-atom distances between body @a index and all others. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			// void calc_pp(unsigned int index);

			/**
			 * @brief Calculate the self-correlation of a body.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_self_correlation(unsigned int index);

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aa(unsigned int n, unsigned int m);

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aw(unsigned int index);

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_ww();

			void combine_self_correlation(unsigned int index);

			void combine_aa(unsigned int n, unsigned int m);

			void combine_aw(unsigned int index);

			void combine_ww();

			/**
			 * @brief Update the compact representation of the coordinates of body @a index.
			 * 
			 * @param index The index of the body to update.
			 */
			void update_compact_representation_body(unsigned int index);

			/**
			 * @brief Update the compact representation of the coordinates of the hydration layer.
			 */
			void update_compact_representation_water();
	};
}