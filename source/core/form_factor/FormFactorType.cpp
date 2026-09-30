// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/FormFactorType.h>

#include <form_factor/FormFactor.h>
#include <form_factor/NeutronFormFactor.h>
#include <settings/ScatteringSettings.h>

using namespace ausaxs::form_factor;

double ausaxs::constants::charge::get_ff_charge(form_factor_t type) {
    // the excluded volume form factor is a normalized shape, so it is shared by both probes
    if (settings::scattering::radiation == settings::scattering::Radiation::Neutron && type != form_factor_t::EXCLUDED_VOLUME) {
        return neutron::protonated::get(type).I0();
    }
    return xray::raw::get(type).I0();
}

double ausaxs::constants::charge::get_ff_charge(form_factor_t type, atom_t fallback_element) {
    // neutron scattering lengths are only tabulated for the form factor types, so OTHER has no better fallback
    if (settings::scattering::radiation == settings::scattering::Radiation::Neutron) {return get_ff_charge(type);}
    if (type == form_factor_t::OTHER) {return ausaxs::constants::charge::nuclear::get_charge(fallback_element);}
    return get_ff_charge(type);
}
