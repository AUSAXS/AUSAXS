#pragma once

#include <utility/Concepts.h>
#include <utility/UtilityFwd.h>
#include <math/MathFwd.h>
#include <io/IOFwd.h>

#include <vector>
#include <utility>

namespace ausaxs::crystal::io {
    struct CrystalReader {
        virtual ~CrystalReader() = default;
        virtual std::pair<Basis3D, std::vector<Vector3<double>>> read(const ausaxs::io::ExistingFile& input) const = 0;
    };
}