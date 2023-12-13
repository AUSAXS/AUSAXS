#pragma once

#include <io/IOFwd.h>
#include <hydrate/Grid.h>

namespace hydrate {
    struct GridReader {
        static grid::Grid read(const io::ExistingFile& filename);
    };
}