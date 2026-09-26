// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/HistFwd.h>
#include <hist/distance_calculator/DistanceCalculatorFwd.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <settings/ExvSettings.h>
#include <utility/observer_ptr.h>

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

/**
 * @brief What the histogram managers build from their distance calculations, shared between the weighted and the form
 *        factor-resolved variants of each manager.
 */
namespace ausaxs::hist::detail {
	/**
	 * @brief The pairwise distance distributions of a histogram manager, before any excluded volume accounting.
	 *        With form_factors, the atomic distributions are resolved by form factor type.
	 */
	template<bool weighted_bins, bool form_factors>
	struct ManagerDistributions {
		using aa_t = std::conditional_t<form_factors,
			typename GenericDistribution3D<weighted_bins, Shape::Triangular>::type, // unordered (ff_type1, ff_type2), distance
			typename GenericDistribution1D<weighted_bins>::type
		>;
		using aw_t = std::conditional_t<form_factors,
			typename GenericDistribution2D<weighted_bins>::type,                     // ff_type, distance
			typename GenericDistribution1D<weighted_bins>::type
		>;
		using ww_t = typename GenericDistribution1D<weighted_bins>::type;

		aa_t p_aa;
		aw_t p_aw;
		ww_t p_ww;
		ww_t p_tot;

		/**
		 * @brief Keep only the first @a bins distance bins of every distribution.
		 */
		void resize(int bins) {
			p_aa.resize(bins);
			p_aw.resize(bins);
			p_ww.resize(bins);
			p_tot.resize(bins);
		}
	};

	/**
	 * @brief Add every histogram of @a from into @a to, which must have the same shape. This works for any dimension,
	 *        since the distance axis is always innermost and contiguous.
	 */
	template<typename T>
	void add_to(T& to, const T& from) {
		std::transform(to.begin(), to.end(), from.begin(), to.begin(), std::plus<>());
	}

	/**
	 * @brief The sum of the results @a ids of @a store, which are all of type @a T.
	 */
	template<typename T, typename Store>
	T sum_results(const Store& store, const std::vector<int>& ids) {
		assert(!ids.empty() && "sum_results: expected at least one result.");
		T total = store.template get<T>(ids.front());
		for (std::size_t i = 1; i < ids.size(); ++i) {add_to(total, store.template get<T>(ids[i]));}
		return total;
	}

	/**
	 * @brief Apply @a op to every bin of @a total and the same bin of each histogram of @a result, i.e. of each class or class
	 *        pair; @a result must span as many bins as @a total. This is how the partial managers keep their total up to date.
	 */
	template<typename Total, typename Result, typename Op>
	void fold_classes(Total& total, const Result& result, Op op) {
		auto bins = static_cast<std::ptrdiff_t>(total.size());
		auto out = total.begin();
		std::ptrdiff_t i = 0;
		for (const auto& bin : result) {
			out[i] = op(out[i], bin);
			if (++i == bins) {i = 0;}
		}
	}

	/**
	 * @brief Export the results @a aa, @a aw and @a ww of @a store, sum them to the total, and trim all four to the bins holding
	 *        anything. The store must have run, and the results must have the shapes of ManagerDistributions.
	 */
	template<bool weighted_bins, bool form_factors>
	ManagerDistributions<weighted_bins, form_factors> export_distributions(distance_calculator::HistogramStore<weighted_bins>& store, int aa, int aw, int ww);

	/**
	 * @brief The composite histogram of the explicit excluded volume models: @a method picks FoXS, Pepsi or CRYSOL, and any other
	 *        method gets the Fraser model.
	 */
	template<bool weighted_bins>
	std::unique_ptr<ICompositeDistanceHistogram> make_explicit_histogram(
		ManagerDistributions<weighted_bins, true>&& distributions, settings::exv::ExvMethod method, observer_ptr<const data::Molecule> protein
	);

	/**
	 * @brief The composite histogram of @a distributions. Weighted ones get the simple model, which needs no @a method.
	 *        Form factor-resolved ones get the average excluded volume model if @a method is ExvMethod::Average, and an explicit
	 *        model otherwise, see make_explicit_histogram.
	 */
	template<bool weighted_bins, bool form_factors>
	std::unique_ptr<ICompositeDistanceHistogram> make_histogram(
		ManagerDistributions<weighted_bins, form_factors>&& distributions, settings::exv::ExvMethod method, observer_ptr<const data::Molecule> protein
	);
}
