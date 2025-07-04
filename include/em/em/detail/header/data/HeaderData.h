// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <unordered_map>

namespace ausaxs::em::detail::header {
    enum class DataType {
        int8,       // int8 --> short int
        int16,      // int16 --> int
        uint8,      // uint8 --> short unsigned int
        uint16,     // uint16 --> unsigned int
        float16,    // float16 --> float
        float32,    // float32 --> float
        complex32,  // complex32 (2x 16bit int)
        complex64,  // complex64 (2x 32bit float)
        NONE        // error type
    };

    extern std::unordered_map<em::detail::header::DataType, unsigned int> byte_sizes;
}