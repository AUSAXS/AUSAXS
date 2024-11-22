#pragma once

#include <crystal/io/CrystalReader.h>

namespace ausaxs::crystal::io {
    struct PDBReader : public CrystalReader {
        ~PDBReader() override = default;
        std::pair<Basis3D, std::vector<Vector3<double>>> read(const ausaxs::io::ExistingFile& input) const override;
    };
}