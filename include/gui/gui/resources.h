// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/ExistingFile.h>

#include <array>
#include <fstream>

using namespace ausaxs;

namespace resources {
    extern const std::array<unsigned char, 24576> resource_file;

    inline io::ExistingFile generate_resource_file() {
        io::File file("resources/elements_basic.ttf");
        if (file.exists()) {return file;}
        else {file.create();}
        std::ofstream out(file.path(), std::ios::binary);
        out.write(reinterpret_cast<const char*>(resource_file.data()), resource_file.size());
        return file;
    }
}