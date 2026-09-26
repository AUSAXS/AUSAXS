// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMT.h>

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool wb>
HistogramManagerMT<wb>::~HistogramManagerMT() = default;

template<bool wb>
std::unique_ptr<DistanceHistogram> HistogramManagerMT<wb>::calculate() {return calculate_all();}

template<bool wb>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMT<wb>::calculate_all() {
    logging::log("HistogramManagerMT::calculate: starting calculation");
    // the simple model has no excluded volume method to choose
    return hist::detail::make_histogram(this->compute_distributions(), settings::exv::ExvMethod::Simple, this->protein);
}

template class hist::HistogramManagerMT<false>;
template class hist::HistogramManagerMT<true>;