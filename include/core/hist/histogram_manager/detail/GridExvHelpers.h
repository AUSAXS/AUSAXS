// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/histogram_manager/detail/ManagerResults.h>

#include <memory>
#include <vector>

namespace ausaxs::hist {
	class ScatteringModelGrid;
}

namespace ausaxs::hist::detail::grid_exv {
	/**
	 * @brief The form factor-resolved atomic distributions, which the grid-based managers extend with their own excluded volume.
	 */
	struct AtomicDistributions {
		Distribution3D<hist::Shape::Triangular> p_aa;
		Distribution2D p_aw;
		Distribution1D p_ww;
		std::vector<double> self_aa;
		double self_ww;
		WeightedDistribution1D p_tot;

		/**
		 * @brief Convert the distributions computed by a form factor manager.
		 */
		static AtomicDistributions from(ManagerDistributions<true, true>&& d);

		/**
		 * @brief Grow every distribution to at least @a bins bins, and return the number of bins they then all have.
		 */
		int grow(int bins);
	};

	/**
	 * @brief Write the grid-based @a p_ax, @a p_wx, and @a p_xx into the excluded volume slot of @a atomic.
	 *        Each of these must span at least as many bins as @a atomic.
	 */
	std::unique_ptr<ScatteringModelGrid> splice(
		AtomicDistributions atomic, const WeightedDistribution2D& p_ax, const WeightedDistribution1D& p_wx, WeightedDistribution1D p_xx
	);
}
