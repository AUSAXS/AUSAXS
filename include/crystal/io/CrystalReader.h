#pragma once

#include <utility/Concepts.h>

#include <vector>
#include <utility>

template<numeric T> class Vector3;
class Basis3D;
namespace io {class ExistingFile;}
namespace crystal::io {
    struct CrystalReader {
        virtual ~CrystalReader() = default;
        virtual std::pair<Basis3D, std::vector<Vector3<double>>> read(const ::io::ExistingFile& input) const = 0;
    };
}