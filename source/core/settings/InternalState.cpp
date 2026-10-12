// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/InternalState.h>

#include <constants/Constants.h>
#include <constants/ConstantsAxes.h>
#include <settings/ScatteringSettings.h>
#include <utility/Exceptions.h>

using namespace ausaxs::settings;

double internal_state::inv_bin_width = 1./constants::axes::d_axis.width();
bool internal_state::prefer_partial_manager = false;
bool internal_state::allow_decorrelate_atom_order = true;

double internal_state::solvent_density() {
    using Radiation = scattering::Radiation;
    switch (scattering::radiation.value) {
        case Radiation::XRay:    return scattering::xray_solvent_density.value;
        case Radiation::Neutron: return scattering::neutron_solvent_density.value;
        default: throw ausaxs::except::unexpected("settings::internal_state::solvent_density: Unknown radiation type.");
    }
}
