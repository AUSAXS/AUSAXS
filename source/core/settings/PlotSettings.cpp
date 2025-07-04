// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/PlotSettings.h>
#include <settings/SettingsIORegistry.h>

using namespace ausaxs;

std::string settings::plots::format = "png";
std::vector<double> settings::plots::contour = {};

namespace ausaxs::settings::io {
    settings::io::SettingSection plot_section("Plot", {
        settings::io::create(plots::format, "format")
    });
}