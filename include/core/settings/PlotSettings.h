// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>
#include <settings/SettingsHelper.h>

#include <string>
#include <vector>

namespace ausaxs::settings {
    struct EXPORT plots {
        static Setting<std::string> format;           // The output format. //! can maybe be merged into general settings
        static Setting<std::vector<double>> contour;  // The contour levels for the image plots. //! unused?
    };
}