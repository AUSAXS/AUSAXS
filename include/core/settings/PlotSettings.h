#pragma once

#include <string>
#include <vector>

namespace settings::plots {
    extern std::string format;           // The output format. //! can maybe be merged into general settings
    extern std::vector<double> contour;  // The contour levels for the image plots. //! unused?
}