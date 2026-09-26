// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>

#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/intensity_calculator/foxs/CompositeDistanceHistogramFoXS.h>
#include <hist/intensity_calculator/pepsi/CompositeDistanceHistogramPepsi.h>
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
    auto raw = this->compute_distributions();

    switch (settings::exv::exv_method) {
        case settings::exv::ExvMethod::FoXS:
            return std::make_unique<CompositeDistanceHistogramFoXS>(
                Distribution3D<hist::Shape::Triangular>(std::move(raw.p_aa)), 
                Distribution2D(std::move(raw.p_aw)), 
                Distribution1D(std::move(raw.p_ww)),
                std::move(raw.p_tot)
            );
        case settings::exv::ExvMethod::Pepsi:
            return std::make_unique<CompositeDistanceHistogramPepsi>(
                Distribution3D<hist::Shape::Triangular>(std::move(raw.p_aa)), 
                Distribution2D(std::move(raw.p_aw)), 
                Distribution1D(std::move(raw.p_ww)),
                std::move(raw.p_tot),
                this->protein
            );
        case settings::exv::ExvMethod::CRYSOL:
            return std::make_unique<CompositeDistanceHistogramCrysol>(
                Distribution3D<hist::Shape::Triangular>(std::move(raw.p_aa)), 
                Distribution2D(std::move(raw.p_aw)), 
                Distribution1D(std::move(raw.p_ww)),
                std::move(raw.p_tot),
                this->protein
            );
        default:
            return std::make_unique<CompositeDistanceHistogramFFExplicit>(
                Distribution3D<hist::Shape::Triangular>(std::move(raw.p_aa)), 
                Distribution2D(std::move(raw.p_aw)), 
                Distribution1D(std::move(raw.p_ww)),
                std::move(raw.p_tot)
            );
    }
}

template class hist::HistogramManagerMTFFExplicit<false>;
template class hist::HistogramManagerMTFFExplicit<true>;
