#pragma once

#include <utility/Container2D.h>

namespace table {
    /**
     * @brief A representation of a table. It uses only a single contiguous vector as storage for better cache locality. 
     */
    using Table = Container2D<double>;
}