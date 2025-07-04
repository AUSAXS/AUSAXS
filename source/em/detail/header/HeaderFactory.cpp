// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <em/detail/header/HeaderFactory.h>
#include <em/detail/header/MapHeader.h>
#include <em/detail/header/MRCHeader.h>
#include <em/detail/header/RECHeader.h>
#include <io/ExistingFile.h>
#include <utility/Exceptions.h>

using namespace ausaxs;

std::unique_ptr<em::detail::header::IMapHeader> em::detail::factory::create_header(const io::ExistingFile& path) {
    if (em::detail::header::MRCHeader::is_mrc(path)) {
        return std::make_unique<em::detail::header::MRCHeader>();
    } else if (em::detail::header::RECHeader::is_rec(path)) {
        return std::make_unique<em::detail::header::RECHeader>();
    } else {
        throw except::io_error("Unsupported map data format: \"" + path.extension() + "\"");
    }
}