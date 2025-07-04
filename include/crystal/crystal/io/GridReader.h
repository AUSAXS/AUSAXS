// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <crystal/io/CrystalReader.h>

namespace ausaxs::crystal::io {
    struct GridReader : public CrystalReader {
        ~GridReader() override = default;
        std::pair<Basis3D, std::vector<Vector3<double>>> read(const ausaxs::io::ExistingFile& input) const override;
    };
}