#pragma once

#include <math/Vector3.h>
#include <utility/Basis3D.h>
#include <io/ExistingFile.h>

namespace crystal::io {
    struct CrystalReader {
        virtual ~CrystalReader() = default;
        virtual std::pair<Basis3D, std::vector<Vector3<double>>> read(const ::io::ExistingFile& input) const = 0;
    };
}