// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/neutron/NeutronFormFactor.h>

#include <form_factor/neutron/NeutronSolvent.h>
#include <settings/ScatteringSettings.h>

using namespace ausaxs;
using namespace ausaxs::form_factor;

neutron::FormFactor neutron::get(form_factor_t type) {
    return get(type, d2o_fraction(settings::scattering::neutron::solvent_density), settings::scattering::neutron::exchange_efficiency);
}
