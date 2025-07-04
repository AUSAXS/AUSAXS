// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>
#include <crystal/io/CrystalIOFwd.h>

#include <memory>

namespace ausaxs::crystal::factory {
    struct CrystalReaderFactory {
        static std::unique_ptr<io::CrystalReader> create(const ausaxs::io::ExistingFile& filename);
    };
}