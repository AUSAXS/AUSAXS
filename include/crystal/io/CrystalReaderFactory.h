#pragma once

#include <memory>

namespace io {class ExistingFile;}
namespace crystal::io {class CrystalReader;}
namespace crystal::factory {
    struct CrystalReaderFactory {
        static std::unique_ptr<io::CrystalReader> create(const ::io::ExistingFile& filename);
    };
}