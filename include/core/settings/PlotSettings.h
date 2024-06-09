#pragma once

#include <string>
#include <vector>

namespace settings {
    namespace plots {
        extern std::string format;           // The output format.
        extern std::vector<double> contour;  // The contour levels for the image plots.
    }
}