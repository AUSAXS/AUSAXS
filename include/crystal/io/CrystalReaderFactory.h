#pragma once

#include <crystal/io/CrystalReader.h>
#include <crystal/io/GridReader.h>
#include <crystal/io/UnitCellReader.h>
#include <io/ExistingFile.h>

namespace crystal::io {
    struct CrystalReaderFactory {
        static std::unique_ptr<CrystalReader> create(const io::ExistingFile& filename);
    };
}