#pragma once

#include <hydrate/Grid.h>

namespace io {class ExistingFile;}
namespace hydrate {
    struct GridReader {
        static grid::Grid read(const io::ExistingFile& filename);
    };
}