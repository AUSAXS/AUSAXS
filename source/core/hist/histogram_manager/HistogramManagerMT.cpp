// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMT.h>

#include <hist/distribution/Distribution1D.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
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
    auto res = this->compute_distributions();
    return std::make_unique<CompositeDistanceHistogram>(
        Distribution1D(std::move(res.p_aa)),
        Distribution1D(std::move(res.p_aw)),
        Distribution1D(std::move(res.p_ww)),
        std::move(res.p_tot)
    );
}

template class hist::HistogramManagerMT<false>;
template class hist::HistogramManagerMT<true>;