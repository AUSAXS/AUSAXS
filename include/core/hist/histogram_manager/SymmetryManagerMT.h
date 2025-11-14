// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/histogram_manager/IHistogramManager.h>
#include <hist/HistFwd.h>
#include <data/DataFwd.h>

#include <memory>

namespace ausaxs::hist {
    template<bool weighted_bins, bool variable_bin_width>
    class SymmetryManagerMT : public IHistogramManager {
        public:
            SymmetryManagerMT(observer_ptr<const data::Molecule> protein);

            std::unique_ptr<hist::DistanceHistogram> calculate() override;

            std::unique_ptr<hist::ICompositeDistanceHistogram> calculate_all() override;

        private:
			observer_ptr<const data::Molecule> protein;

            template<bool contains_waters>
            std::unique_ptr<hist::ICompositeDistanceHistogram> calculate();
    };
}