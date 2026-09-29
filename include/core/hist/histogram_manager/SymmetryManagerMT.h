// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/HistFwd.h>
#include <hist/histogram_manager/IHistogramManager.h>

#include <memory>

namespace ausaxs::hist {
    /**
     * @brief The multithreaded histogram manager for molecules with symmetries, which calculates the whole
     *        histogram in one go. Each symmetric copy is only evaluated once, and the histogram scaled by how often it occurs.
     *
     * @tparam form_factors Whether the atoms are resolved by form factor, see HistogramManagerMTBase.
     */
    template<bool weighted_bins, bool form_factors>
    class SymmetryManagerMTBase : public IHistogramManager {
        public:
            explicit SymmetryManagerMTBase(observer_ptr<const data::Molecule> protein);

            std::unique_ptr<hist::DistanceHistogram> calculate() override;

            std::unique_ptr<hist::ICompositeDistanceHistogram> calculate_all() override;

        private:
			observer_ptr<const data::Molecule> protein;

            template<bool contains_waters>
            std::unique_ptr<hist::ICompositeDistanceHistogram> calculate();
    };

    /**
     * @brief The symmetry manager for the simple excluded volume model, where every atom carries its own weight.
     */
    template<bool weighted_bins>
    using SymmetryManagerMT = SymmetryManagerMTBase<weighted_bins, false>;

    /**
     * @brief The symmetry manager for the form factor-resolved excluded volume models.
     */
    template<bool weighted_bins>
    using SymmetryManagerMTFF = SymmetryManagerMTBase<weighted_bins, true>;
}