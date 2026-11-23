// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/IScatteringModel.h>

#include <hist/distribution/Distribution1D.h>
#include <hist/scattering_model/AmplitudeFormFactors.h>

using namespace ausaxs;

bool hist::IScatteringModel::is_highly_ordered() const {
    return DistanceHistogram::is_highly_ordered(get_aa_counts().get_content());
}

std::vector<double> hist::IScatteringModel::amplitude_transform(
    const AmplitudeScatterers& scatterers, const std::vector<double>& q, const std::optional<Vector3<double>>& centre
) const {
    return hist::amplitude_transform(amplitude_form_factors.get(), scatterers, q, get_parameters(), centre);
}

void hist::IScatteringModel::set_amplitude_form_factors(std::shared_ptr<const AmplitudeFormFactors> form_factors) {
    amplitude_form_factors = std::move(form_factors);
}

hist::ScatteringParameters hist::IScatteringModel::get_parameters() const {return {};}
