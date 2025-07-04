// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <em/detail/header/data/HeaderData.h>

#include <cstdint>

using namespace ausaxs;

namespace ausaxs::em::detail::header {
    std::unordered_map<em::detail::header::DataType, unsigned int> byte_sizes = {
        {em::detail::header::DataType::int8,    sizeof(int8_t)},
        {em::detail::header::DataType::int16,   sizeof(int16_t)},
        {em::detail::header::DataType::uint8,   sizeof(uint8_t)},
        {em::detail::header::DataType::uint16,  sizeof(uint16_t)},
        {em::detail::header::DataType::float16, 2u},
        {em::detail::header::DataType::float32, sizeof(float)},
    };
}