#pragma once

#include <utility/SmartOption.h>

namespace settings {
    namespace plots {
        extern settings::detail::SmartOption<std::string> format;           // The output format.
        extern settings::detail::SmartOption<std::vector<double>> contour;  // The contour levels for the image plots.
    }
}