// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/FormFactorType.h>

#include <form_factor/FormFactor.h>
#include <form_factor/NeutronFormFactor.h>
#include <settings/ScatteringSettings.h>

using namespace ausaxs::form_factor;
using Radiation = ausaxs::settings::scattering::Radiation;

template<Radiation radiation>
double ausaxs::constants::charge::get_ff_charge(form_factor_t type) {
    // the excluded volume form factor is a normalized shape, so it is shared by both probes
    if constexpr (radiation == Radiation::Neutron) {
        if (type != form_factor_t::EXCLUDED_VOLUME) {return neutron::protonated::get(type).I0();}
    }
    return xray::raw::get(type).I0();
}

template<Radiation radiation>
double ausaxs::constants::charge::get_ff_charge(form_factor_t type, atom_t fallback_element) {
    // neutron scattering lengths are only tabulated for the form factor types, so OTHER has no better fallback
    if constexpr (radiation == Radiation::XRay) {
        if (type == form_factor_t::OTHER) {return nuclear::get_charge(fallback_element);}
    }
    return get_ff_charge<radiation>(type);
}

template double ausaxs::constants::charge::get_ff_charge<Radiation::XRay>(form_factor_t);
template double ausaxs::constants::charge::get_ff_charge<Radiation::Neutron>(form_factor_t);
template double ausaxs::constants::charge::get_ff_charge<Radiation::XRay>(form_factor_t, atom_t);
template double ausaxs::constants::charge::get_ff_charge<Radiation::Neutron>(form_factor_t, atom_t);

double ausaxs::constants::charge::get_ff_charge(form_factor_t type) {
    switch (settings::scattering::radiation.value) {
        case Radiation::XRay:    return get_ff_charge<Radiation::XRay>(type);
        case Radiation::Neutron: return get_ff_charge<Radiation::Neutron>(type);
    }
    throw except::unexpected("constants::charge::get_ff_charge: Unknown radiation type.");
}

double ausaxs::constants::charge::get_ff_charge(form_factor_t type, atom_t fallback_element) {
    switch (settings::scattering::radiation.value) {
        case Radiation::XRay:    return get_ff_charge<Radiation::XRay>(type, fallback_element);
        case Radiation::Neutron: return get_ff_charge<Radiation::Neutron>(type, fallback_element);
    }
    throw except::unexpected("constants::charge::get_ff_charge: Unknown radiation type.");
}
