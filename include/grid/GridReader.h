#pragma once

#include <io/IOFwd.h>
#include <grid/Grid.h>

namespace grid {
    struct GridReader {
        static grid::Grid read(const io::ExistingFile& filename);
    };
}