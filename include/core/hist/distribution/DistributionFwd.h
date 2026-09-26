// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/indexers/Shape.h>

namespace ausaxs::hist {
    using utility::indexer::Shape;

    class Distribution1D;
    class Distribution2D;
    template<Shape S = Shape::Square> class Distribution3D;
    class WeightedDistribution1D;
    class WeightedDistribution2D;
    template<Shape S = Shape::Square> class WeightedDistribution3D;
}