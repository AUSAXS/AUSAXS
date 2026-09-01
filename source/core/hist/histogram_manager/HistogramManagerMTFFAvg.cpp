// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <data/Molecule.h>
#include <constants/Constants.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool wb, bool vbw>
HistogramManagerMTFFAvg<wb, vbw>::~HistogramManagerMTFFAvg() = default;

template<bool wb, bool vbw>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFAvg<wb, vbw>::calculate() {return calculate_all();}

template<bool wb, bool vbw>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFAvg<wb, vbw>::calculate_all() {
    logging::log("HistogramManagerMTFFAvg::calculate: starting calculation");
    auto raw = this->compute_raw_distributions();

    // The average excluded volume charge displaced by a single atom.
    double Z_exv_avg = 
        this->protein->size_atom() == 0 
        ? 0 
        : this->protein->get_volume_grid()*constants::charge::density::water/this->protein->size_atom()
    ;

    return std::make_unique<CompositeDistanceHistogramFFAvg>(
        Distribution3D(std::move(raw.p_aa)), 
        Distribution2D(std::move(raw.p_aw)), 
        Distribution1D(std::move(raw.p_ww)), 
        std::move(raw.p_tot),
        Z_exv_avg
    );
}

template class hist::HistogramManagerMTFFAvg<false, false>;
template class hist::HistogramManagerMTFFAvg<false, true>;
template class hist::HistogramManagerMTFFAvg<true, false>;
template class hist::HistogramManagerMTFFAvg<true, true>;
