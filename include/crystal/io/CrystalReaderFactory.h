#pragma once

#include <crystal/io/CrystalReader.h>
#include <crystal/io/GridReader.h>
#include <crystal/io/UnitCellReader.h>

namespace crystal::io {
    struct CrystalReaderFactory {
        static std::unique_ptr<CrystalReader> create(const std::string& filename);
    };
}