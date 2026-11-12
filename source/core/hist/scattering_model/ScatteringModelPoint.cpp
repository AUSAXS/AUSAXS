// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/ScatteringModelPoint.h>

#include <utility>

using namespace ausaxs;
using namespace ausaxs::hist;

ScatteringModelPoint::ScatteringModelPoint(DistanceHistogram&& hist) noexcept : DistanceHistogram(std::move(hist)) {}
