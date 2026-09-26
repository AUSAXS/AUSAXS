// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/WeightedDistribution2D.h>

#include <memory>

namespace ausaxs::hist {
	class CompositeDistanceHistogramFFAvg;
	class CompositeDistanceHistogramFFGrid;
}

namespace ausaxs::hist::detail::grid_exv {
	/**
	 * @brief The distributions of a form factor-averaged result, which the grid-based managers extend with their own excluded volume.
	 */
	struct AtomicDistributions {
		Distribution3D<hist::Shape::Triangular> p_aa;
		Distribution2D p_aw;
		Distribution1D p_ww;
		WeightedDistribution1D p_tot;

		/**
		 * @brief Move the distributions out of @a base, which must not be used afterwards.
		 */
		static AtomicDistributions take(CompositeDistanceHistogramFFAvg& base);

		/**
		 * @brief Grow every distribution to at least @a bins bins, and return the number of bins they then all have.
		 */
		int grow(int bins);
	};

	/**
	 * @brief Replace the excluded volume of the averaged model in @a atomic with the grid-based @a p_ax, @a p_wx, and @a p_xx.
	 *        Each of these must span at least as many bins as @a atomic.
	 */
	std::unique_ptr<CompositeDistanceHistogramFFGrid> splice(
		AtomicDistributions atomic, const WeightedDistribution2D& p_ax, const WeightedDistribution1D& p_wx, WeightedDistribution1D p_xx
	);
}
