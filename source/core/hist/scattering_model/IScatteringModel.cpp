// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/IScatteringModel.h>

#include <hist/distribution/Distribution1D.h>

using namespace ausaxs;

bool hist::IScatteringModel::is_highly_ordered() const {
    return DistanceHistogram::is_highly_ordered(get_aa_counts().get_content());
}
