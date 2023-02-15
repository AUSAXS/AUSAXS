#pragma once

#include <dataset/detail/DatasetConstructor.h>

namespace detail {
    /**
     * @brief Contructs a dataset from a DAT file.
     *        Supports the formats x | y | yerr and x | y | yerr | xerr.
     */
    struct DATConstructor : DatasetConstructor {
        /**
         * @brief Construct a dataset from a DAT file.
         */
        std::shared_ptr<Dataset> construct(std::string path);
    };
}