// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <array>
#include <string>
#include <string_view>

namespace ausaxs::detail::qunit {
    /**
     * @brief Check if a header line already names a q-unit.
     */
    inline bool header_names_a_unit(std::string_view header) {
        constexpr std::array<std::string_view, 4> units = {"[A]", "[A^-1]", "[nm]", "[nm^-1]"};
        for (std::size_t start = 0; start < header.size();) {
            std::size_t end = header.find_first_of(" \t", start);
            std::string_view token = header.substr(start, end == std::string_view::npos ? std::string_view::npos : end-start);
            for (const auto& unit : units) {
                if (token == unit) {return true;}
            }
            if (end == std::string_view::npos) {break;}
            start = end+1;
        }
        return false;
    }

    /**
     * @brief Get the unit annotation line for a file written by one of our own save methods.
     *
     * All q-values are stored in inverse Ångström internally, so that is what we always write. Without the annotation the reader falls 
     * back to guessing the unit from the magnitude of the q-values, which silently rescales our own output by a factor 10 whenever the 
     * nanometre default is configured.
     */
    inline std::string unit_line(std::string_view header) {
        return header_names_a_unit(header) ? "" : "[A]\n";
    }
}
