#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <settings/HistogramSettings.h>

#include "hist/intensity_calculator/DistanceHistogram.h"
#include "hist/hist_test_helper.h"

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::data;

TEST_CASE("Custom bin width: Managers use correct bin width") {
    settings::axes::bin_width = GENERATE(0.1, 0.5, 1.0);
}