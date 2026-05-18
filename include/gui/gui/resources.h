// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/ExistingFile.h>

#include <array>
#include <fstream>

using namespace ausaxs;

namespace resources {
    /// @brief The embedded font file (elements_basic.ttf), stored as a raw byte array.
    extern const std::array<unsigned char, 24576> resource_file;

    /**
     * @brief Write the embedded font to disk and return a handle to it.
     *
     * If the resource file already exists on disk it is reused as-is; otherwise it is created from
     * the embedded @ref resource_file byte array.
     */
    inline io::ExistingFile generate_resource_file() {
        io::File file("resources/elements_basic.ttf");
        if (file.exists()) {return file;}
        else {file.create();}
        std::ofstream out(file.path(), std::ios::binary);
        out.write(reinterpret_cast<const char*>(resource_file.data()), resource_file.size());
        return file;
    }
}