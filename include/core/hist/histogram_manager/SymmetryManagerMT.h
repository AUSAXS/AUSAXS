// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/HistFwd.h>
#include <hist/histogram_manager/IHistogramManager.h>
#include <settings/ExvSettings.h>

#include <memory>

namespace ausaxs::hist {
    /**
     * @brief Common machinery of the multithreaded histogram managers for molecules with symmetries, which calculate the whole
     *        histogram in one go. Each symmetric copy is only evaluated once, and the histogram scaled by how often it occurs.
     *
     * @tparam form_factors Whether the atoms are resolved by form factor, see HistogramManagerMTBase.
     */
    template<bool weighted_bins, bool form_factors>
    class SymmetryManagerMTBase : public IHistogramManager {
        public:
            /**
             * @param exv_method The excluded volume model the form factor-resolved result is built for; see detail::make_histogram.
             */
            SymmetryManagerMTBase(observer_ptr<const data::Molecule> protein, settings::exv::ExvMethod exv_method);

            std::unique_ptr<hist::DistanceHistogram> calculate() override;

            std::unique_ptr<hist::ICompositeDistanceHistogram> calculate_all() override;

        private:
			observer_ptr<const data::Molecule> protein;
            settings::exv::ExvMethod exv_method;

            template<bool contains_waters>
            std::unique_ptr<hist::ICompositeDistanceHistogram> calculate();
    };

    /**
     * @brief The symmetry manager for the simple excluded volume model, where every atom carries its own weight.
     */
    template<bool weighted_bins>
    // NOLINTNEXTLINE - the destructor is virtual through the dependent base, which the check cannot see on the template pattern
    class SymmetryManagerMT : public SymmetryManagerMTBase<weighted_bins, false> {
        public:
            explicit SymmetryManagerMT(observer_ptr<const data::Molecule> protein)
                : SymmetryManagerMTBase<weighted_bins, false>(protein, settings::exv::ExvMethod::Simple) {}
    };

    /**
     * @brief The symmetry manager for the form factor-resolved excluded volume models.
     */
    template<bool weighted_bins>
    // NOLINTNEXTLINE - the destructor is virtual through the dependent base, which the check cannot see on the template pattern
    class SymmetryManagerMTFF : public SymmetryManagerMTBase<weighted_bins, true> {
        public:
            explicit SymmetryManagerMTFF(observer_ptr<const data::Molecule> protein, settings::exv::ExvMethod exv_method = settings::exv::exv_method)
                : SymmetryManagerMTBase<weighted_bins, true>(protein, exv_method) {}
    };
}