#pragma once

namespace BS {class thread_pool;}

#include <hist/PartialHistogramManager.h>

#include <memory>

namespace hist {
	/**
	 * @brief A multi-threaded smart distance calculator.
	 */
	class PartialHistogramManagerMT : public PartialHistogramManager {
		public:
			PartialHistogramManagerMT(Protein* protein);

			PartialHistogramManagerMT(PartialHistogramManager&);

			~PartialHistogramManagerMT() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			Histogram calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			ScatteringHistogram calculate_all() override;

		private:
			std::unique_ptr<BS::thread_pool> pool;

			/**
			 * @brief Initialize this object. The internal distances between atoms in each body is constant and cannot change. 
			 *        They are unaffected by both rotations and translations, and so we precalculate them. 
			 */
			void initialize();

			/**
			 * @brief Calculate the atom-atom distances between body @a index and all others. 
			 */
			void calc_pp(unsigned int index);

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 */
			void calc_pp(unsigned int n, unsigned int m);

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 */
			void calc_hp(unsigned int index);

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 */
			void calc_hh();

			/**
			 * @brief Calculate the self-correlation of a body. 
			 */
			void calc_self_correlation(unsigned int index);
	};
}