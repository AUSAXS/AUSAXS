#pragma once

#include <string>
#include <vector>

namespace ausaxs::settings {
    struct plots {
        static std::string format;           // The output format. //! can maybe be merged into general settings
        static std::vector<double> contour;  // The contour levels for the image plots. //! unused?
    };
}