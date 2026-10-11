// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/PlotSettings.h>

#include <settings/SettingsIORegistry.h>

using namespace ausaxs;

settings::Setting<std::string> settings::plots::format{"png"};
settings::Setting<std::vector<double>> settings::plots::contour{{}};

namespace {
    using namespace ausaxs::settings;
    settings::io::SettingSection plot_section("Plot", {
        settings::io::create(plots::format, "format")
    });
}