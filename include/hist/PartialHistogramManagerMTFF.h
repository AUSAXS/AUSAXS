#pragma once

#include <hist/HistogramManager.h>
#include <hist/detail/MasterHistogram.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <utility/Container1D.h>
#include <utility/Container2D.h>
#include <utility/Container3D.h>

#include <memory>

namespace BS {
	class thread_pool;
	template <typename T> class multi_future;
}
namespace hist {
	/**
	 * @brief A multi-threaded smart distance calculator.
	 */
	class PartialHistogramManagerMTFF : public HistogramManager {
		public:
			PartialHistogramManagerMTFF(Protein* protein);

			~PartialHistogramManagerMTFF() override;

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
			std::mutex master_hist_mutex;

			detail::MasterHistogram master;                       	// the current total histogram
			std::vector<detail::CompactCoordinatesFF> coords_p;
			detail::CompactCoordinatesFF coords_h;
			Container3D<detail::PartialHistogram> partials_pp;		// x: form factor index, y: body index, z: body index 
			Container2D<detail::HydrationHistogram> partials_hp; 	// x: form factor index, y: body index
			Container1D<detail::HydrationHistogram> partials_hh;	// x: form factor index

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
			BS::multi_future<std::vector<double>> calc_self_correlation(unsigned int index);

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			BS::multi_future<std::vector<double>> calc_pp(unsigned int n, unsigned int m);

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			BS::multi_future<std::vector<double>> calc_hp(unsigned int index);

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 * 		  This only adds jobs to the thread pool, and does not wait for them to complete.
			 */
			BS::multi_future<std::vector<double>> calc_hh();

			void combine_self_correlation(unsigned int index, BS::multi_future<std::vector<double>>& futures);

			void combine_pp(unsigned int n, unsigned int m, BS::multi_future<std::vector<double>>& futures);

			void combine_hp(unsigned int index, BS::multi_future<std::vector<double>>& futures);

			void combine_hh(BS::multi_future<std::vector<double>>& futures);

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