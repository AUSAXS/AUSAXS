#pragma once

#include <io/IOFwd.h>
#include <em/detail/header/HeaderFwd.h>

#include <memory>

namespace em::detail::factory {
    std::unique_ptr<em::detail::header::IMapHeader> create_header(const io::ExistingFile& path);
}