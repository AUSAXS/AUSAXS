#pragma once

#include <crystal/io/CrystalReader.h>
#include <crystal/io/GridReader.h>
#include <crystal/io/UnitCellReader.h>

namespace crystal::factory {
    struct CrystalReaderFactory {
        static std::unique_ptr<io::CrystalReader> create(const ::io::ExistingFile& filename);
    };
}