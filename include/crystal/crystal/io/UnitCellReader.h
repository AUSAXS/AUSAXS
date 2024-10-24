#pragma once

#include <crystal/io/CrystalReader.h>

namespace ausaxs::crystal::io {
    /**
     * @brief Read a unit cell from a file. 
     * 
     * The file format specification is as follows:
     *      1. BASIS
     *      2. The x basis vector
     *      3. The y basis vector
     *      4. The z basis vector
     *      5. (empty)
     *      6. CRYSTALDATA
     *      7. Start of data section. Format: x y z
     *      8. End of data section.
     */
    struct UnitCellReader : public CrystalReader {
        ~UnitCellReader() override = default;
        std::pair<Basis3D, std::vector<Vector3<double>>> read(const ::io::ExistingFile& input) const override;
    };
}