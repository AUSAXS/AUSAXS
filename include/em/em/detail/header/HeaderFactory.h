#pragma once

#include <io/IOFwd.h>

#include <memory>

namespace em::detail::header {class MapHeader;}
namespace em::detail::factory {
    std::unique_ptr<em::detail::header::MapHeader> create_header(const io::ExistingFile& path);
}