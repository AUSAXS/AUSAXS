#pragma once

#include <dataset/detail/DatasetConstructor.h>

namespace detail {
    /**
     * @brief Contructs a dataset from a GROMACS XVG file.
     */
    struct XVGConstructor : DatasetConstructor {
        /**
         * @brief Construct a dataset from a GROMACS XVG file.
         */
        std::shared_ptr<Dataset> construct(std::string path);
    };
}