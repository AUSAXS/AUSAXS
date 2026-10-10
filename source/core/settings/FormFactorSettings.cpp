// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/FormFactorSettings.h>

#include <settings/SettingsIORegistry.h>

using namespace ausaxs;

int settings::form_factor::max_types = 15;
double settings::form_factor::min_fraction = 0.001;

namespace {
    using namespace ausaxs::settings;
    settings::io::SettingSection form_factor_section("Form factors", {
        settings::io::create(settings::form_factor::max_types, "max_form_factor_types"),
        settings::io::create(settings::form_factor::min_fraction, "min_form_factor_fraction")
    });
}
