// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/IHistogramManager.h>
#include <hist/scattering_model/ScatteringAmplitude.h>
#include <utility/Exceptions.h>

using namespace ausaxs;

std::unique_ptr<hist::ScatteringAmplitude> hist::IHistogramManager::get_scattering_amplitude() const {
    throw except::invalid_argument(
        "IHistogramManager::get_scattering_amplitude: The scattering amplitude is not available for this histogram manager. "
        "Molecules with symmetries are not supported."
    );
}
