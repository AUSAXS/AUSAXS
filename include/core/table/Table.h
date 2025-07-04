// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/Container2D.h>

namespace ausaxs::table {
    /**
     * @brief A representation of a table. It uses only a single contiguous vector as storage for better cache locality. 
     */
    using Table = container::Container2D<double>;
}