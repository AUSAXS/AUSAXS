#pragma once

#include <crystal/io/CrystalReader.h>

namespace crystal::io {
    struct GridReader : public CrystalReader {
        ~GridReader() override = default;
        std::pair<Basis3D, std::vector<Vector3<double>>> read(const ::io::ExistingFile& input) const override;
    };
}