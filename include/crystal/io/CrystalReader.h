#pragma once

#include <math/Vector3.h>
#include <utility/Basis3D.h>

namespace crystal::io {
    struct CrystalReader {
        virtual ~CrystalReader() = default;
        virtual std::pair<Basis3D, std::vector<Vector3<double>>> read(const std::string& input) const = 0;
    };
}