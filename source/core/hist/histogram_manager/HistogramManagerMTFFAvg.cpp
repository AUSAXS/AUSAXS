// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool wb>
HistogramManagerMTFFAvg<wb>::~HistogramManagerMTFFAvg() = default;

template<bool wb>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFAvg<wb>::calculate() {return calculate_all();}

template<bool wb>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFAvg<wb>::calculate_all() {
    logging::log("HistogramManagerMTFFAvg::calculate: starting calculation");
    return hist::detail::make_histogram(this->compute_distributions(), settings::exv::ExvMethod::Average, this->protein);
}

template class hist::HistogramManagerMTFFAvg<false>;
template class hist::HistogramManagerMTFFAvg<true>;
