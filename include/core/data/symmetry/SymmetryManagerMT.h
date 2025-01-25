#pragma once

#include <data/DataFwd.h>
#include <hist/HistFwd.h>

#include <memory>

namespace ausaxs::symmetry {
    class SymmetryManagerMT {
        public:
            template<bool use_weighted_distribution>
            std::unique_ptr<hist::ICompositeDistanceHistogram> calculate(const data::Molecule& protein);
        
        private:
            template<bool use_weighted_distribution, bool contains_waters>
            std::unique_ptr<hist::ICompositeDistanceHistogram> calculate(const data::Molecule& protein);
    };
}