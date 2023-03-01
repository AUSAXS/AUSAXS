#pragma once

#include <hydrate/Grid.h>

namespace hydrate {
    struct GridReader {
        static Grid read(const std::string& filename);
    };
}