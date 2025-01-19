#pragma once

#include <hist/histogram_manager/PartialHistogramManagerMT.h>

namespace ausaxs::hist {
	/**
	 * @brief A multi-threaded smart distance calculator which efficiently calculates the simple distance histogram. 
	 */
    template<bool use_weighted_distribution> 
	class PartialSymmetryManagerMT : public PartialHistogramManagerMT<use_weighted_distribution> {
		using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
		public:
			using PartialHistogramManagerMT<use_weighted_distribution>::PartialHistogramManagerMT;
			virtual ~PartialSymmetryManagerMT() override;

			/**
			 * @brief Initialize this object. The internal distances between atoms in each body is constant and cannot change. 
			 *        They are unaffected by both rotations and translations, and so we precalculate them. 
			 */
			void initialize() override;

			/**
			 * @brief Calculate the atom-atom distances between body @a index and all others. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			// void calc_pp(unsigned int index);

			/**
			 * @brief Calculate the self-correlation of a body.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_self_correlation(unsigned int index) override;

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aa(unsigned int n, unsigned int m) override;

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_aw(unsigned int index) override;

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			void calc_ww() override;

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