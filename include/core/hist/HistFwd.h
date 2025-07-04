// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::hist {
    class ICompositeDistanceHistogram;
    class ICompositeDistanceHistogramExv;
    class IHistogramManager;
    class DistanceHistogram;
    class Histogram;
    class Histogram2D;

    /**
     * @brief A ScatteringProfile is just a (q, I(q)) histogram. 
     */    
    using ScatteringProfile = Histogram;
}