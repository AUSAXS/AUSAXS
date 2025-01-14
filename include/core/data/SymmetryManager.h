#pragma once

#include <data/DataFwd.h>
#include <hist/HistFwd.h>

#include <memory>

namespace ausaxs::hist::detail {
    class SymmetryManager {
        public:
            template<bool use_weighted_distribution>
            std::unique_ptr<ICompositeDistanceHistogram> calculate(const data::Molecule& protein);
        
        private:
            template<bool use_weighted_distribution, bool contains_waters>
            std::unique_ptr<ICompositeDistanceHistogram> calculate(const data::Molecule& protein);
    };
}