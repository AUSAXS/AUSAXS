// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>

#include <hist/scattering_model/IScatteringModel.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool wb>
HistogramManagerMTFFExplicit<wb>::~HistogramManagerMTFFExplicit() = default;

template<bool wb>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFExplicit<wb>::calculate() {return calculate_all();}

template<bool wb>
std::unique_ptr<IScatteringModel> HistogramManagerMTFFExplicit<wb>::calculate_all() {
    logging::log("HistogramManagerMTFFExplicit::calculate: starting calculation");
    return hist::detail::make_explicit_model(this->compute_distributions(), this->protein);
}

template class hist::HistogramManagerMTFFExplicit<false>;
template class hist::HistogramManagerMTFFExplicit<true>;
