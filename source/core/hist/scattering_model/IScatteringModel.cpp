// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/IScatteringModel.h>

#include <hist/distribution/Distribution1D.h>
#include <hist/scattering_model/AmplitudeHistogram.h>
#include <math/Vector3.h>
#include <utility/Exceptions.h>

using namespace ausaxs;

bool hist::IScatteringModel::is_highly_ordered() const {
    return DistanceHistogram::is_highly_ordered(get_aa_counts().get_content());
}

std::vector<double> hist::IScatteringModel::amplitude_transform(
    const AmplitudeScatterers& scatterers, const std::vector<double>& q, const std::optional<Vector3<double>>& centre
) const {
    auto g = [this, form_factors = scatterers.form_factors] (double q) {return amplitude_form_factors(form_factors, q);};
    auto origin = centre.has_value() ? *centre : AmplitudeHistogram::scattering_centre(scatterers, g(0));
    return AmplitudeHistogram(scatterers, origin).transform(q, g);
}

std::vector<double> hist::IScatteringModel::amplitude_form_factors(bool /*form_factors*/, double /*q*/) const {
    throw except::invalid_argument(
        "IScatteringModel::amplitude_transform: The scattering amplitude is only available for the excluded volume models "
        "\"simple\", \"none\", \"vacuum\" and \"fraser\"."
    );
}
