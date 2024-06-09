#pragma once

#include <io/IOFwd.h>
#include <crystal/io/CrystalIOFwd.h>

#include <memory>

namespace crystal::factory {
    struct CrystalReaderFactory {
        static std::unique_ptr<io::CrystalReader> create(const ::io::ExistingFile& filename);
    };
}