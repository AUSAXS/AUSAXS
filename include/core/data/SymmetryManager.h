#pragma once

#include <math/Matrix.h>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>

namespace ausaxs::hist::detail {
    class SymmetryManager {
        public:
            template<bool use_weighted_distribution>
            std::unique_ptr<CompositeDistanceHistogram> calculate(const data::Molecule& protein);
    };
}