// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/foxs/CompositeDistanceHistogramFoXS.h>
#include <hist/intensity_calculator/pepsi/CompositeDistanceHistogramPepsi.h>
#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <data/Molecule.h>
#include <settings/ExvSettings.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool wb, bool vbw>
HistogramManagerMTFFExplicit<wb, vbw>::~HistogramManagerMTFFExplicit() = default;

template<bool wb, bool vbw>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFExplicit<wb, vbw>::calculate() {return calculate_all();}

template<bool wb, bool vbw>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFExplicit<wb, vbw>::calculate_all() {
    logging::log("HistogramManagerMTFFExplicit::calculate: starting calculation");
    auto raw = this->compute_raw_distributions();

    // No excluded volume distance data is synthesized; each of the calculators below bakes a
    // per-species excluded volume form factor into its own form factor lookup tables instead.
    switch (settings::exv::exv_method) {
        case settings::exv::ExvMethod::FoXS:
            return std::make_unique<CompositeDistanceHistogramFoXS>(
                Distribution3D(std::move(raw.p_aa)), 
                Distribution2D(std::move(raw.p_aw)), 
                Distribution1D(std::move(raw.p_ww)),
                std::move(raw.p_tot)
            );
        case settings::exv::ExvMethod::Pepsi:
            return std::make_unique<CompositeDistanceHistogramPepsi>(
                Distribution3D(std::move(raw.p_aa)), 
                Distribution2D(std::move(raw.p_aw)), 
                Distribution1D(std::move(raw.p_ww)),
                std::move(raw.p_tot),
                form_factor::ExvTableManager::get_average_displaced_volume(this->protein)
            );
        case settings::exv::ExvMethod::CRYSOL:
            return std::make_unique<CompositeDistanceHistogramCrysol>(
                Distribution3D(std::move(raw.p_aa)), 
                Distribution2D(std::move(raw.p_aw)), 
                Distribution1D(std::move(raw.p_ww)),
                std::move(raw.p_tot),
                form_factor::ExvTableManager::get_average_displaced_volume(this->protein)
            );
        default:
            return std::make_unique<CompositeDistanceHistogramFFExplicit>(
                Distribution3D(std::move(raw.p_aa)), 
                Distribution2D(std::move(raw.p_aw)), 
                Distribution1D(std::move(raw.p_ww)),
                std::move(raw.p_tot)
            );
    }
}

template class hist::HistogramManagerMTFFExplicit<false, false>;
template class hist::HistogramManagerMTFFExplicit<false, true>;
template class hist::HistogramManagerMTFFExplicit<true, false>;
template class hist::HistogramManagerMTFFExplicit<true, true>;
