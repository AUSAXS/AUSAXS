// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <settings/ExvSettings.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool wb>
HistogramManagerMTFFExplicit<wb>::~HistogramManagerMTFFExplicit() = default;

template<bool wb>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFExplicit<wb>::calculate() {return calculate_all();}

template<bool wb>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFExplicit<wb>::calculate_all() {
    logging::log("HistogramManagerMTFFExplicit::calculate: starting calculation");
    return hist::detail::make_explicit_histogram(this->compute_distributions(), settings::exv::exv_method, this->protein);
}

template class hist::HistogramManagerMTFFExplicit<false>;
template class hist::HistogramManagerMTFFExplicit<true>;
