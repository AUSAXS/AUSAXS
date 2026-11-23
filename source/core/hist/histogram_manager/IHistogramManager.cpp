// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/IHistogramManager.h>
#include <hist/scattering_model/AmplitudeHistogram.h>
#include <utility/Exceptions.h>

using namespace ausaxs;

hist::AmplitudeScatterers hist::IHistogramManager::get_amplitude_scatterers() const {
    throw except::invalid_argument(
        "IHistogramManager::get_amplitude_scatterers: The scattering amplitude is not available for this histogram manager. "
        "Molecules with symmetries are not supported."
    );
}
