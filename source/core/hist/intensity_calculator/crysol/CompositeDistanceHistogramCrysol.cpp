// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>

#include <form_factor/lookup/ExvTableManager.h>
#include <settings/ExvSettings.h>

#include <cmath>
#include <numbers>

using namespace ausaxs;
using namespace ausaxs::hist;

namespace {
    /**
     * @brief Switch to the Traube volumes used by CRYSOL, and get the average displaced volume per atom of @a molecule.
     *        The switch must happen first, since the average volume depends on the volume set.
     */
    double use_traube_volumes(observer_ptr<const data::Molecule> molecule) {
        // only assign when needed, since every assignment rebuilds the form factor tables
        if (settings::exv::exv_set != settings::exv::ExvSet::Traube) {settings::exv::exv_set = settings::exv::ExvSet::Traube;}
        return form_factor::ExvTableManager::get_average_displaced_volume(molecule);
    }
}

CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot,
    observer_ptr<const data::Molecule> molecule
) : CompositeDistanceHistogramFFExplicit(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), average_displaced_V(use_traube_volumes(molecule)) {}

CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot,
    observer_ptr<const data::Molecule> molecule
) : CompositeDistanceHistogramFFExplicit(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), average_displaced_V(use_traube_volumes(molecule)) {}

double CompositeDistanceHistogramCrysol::exv_factor(double q, double cx, double avg_displaced_V) {
    // G(q) factor from CRYSOL: https://doi.org/10.1107/S0021889895007047
    double c = constexpr_math::pow(avg_displaced_V, 2./3)/(4*std::numbers::pi);
    return std::pow(cx, 3)*std::exp(-c*(std::pow(cx, 2) - 1)*q*q);
}

double CompositeDistanceHistogramCrysol::exv_factor(double q) const {
    return exv_factor(q, free_params.cx, average_displaced_V);
}

Limit CompositeDistanceHistogramCrysol::get_excluded_volume_scaling_factor_limits() const {
    return {0.8, 1.265};
}
