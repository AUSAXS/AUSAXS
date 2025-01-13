#pragma once

#include <data/DataFwd.h>
#include <hist/HistFwd.h>

#include <memory>

namespace ausaxs::hist::detail {
    class SymmetryManager {
        public:
            template<bool use_weighted_distribution>
            std::unique_ptr<ICompositeDistanceHistogram> calculate(const data::Molecule& protein);
    };
}