#pragma once

#include <crystal/io/CrystalReader.h>

namespace crystal::io {
    struct PDBReader : public CrystalReader {
        ~PDBReader() override = default;
        std::pair<Basis3D, std::vector<Vector3<double>>> read(const std::string& input) const override;
    };
}