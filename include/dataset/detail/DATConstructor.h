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
         * 
         * @param path The path to the file.
         * @param expected_cols The expected number of columns. Any additional columns will be ignored. If 0, all columns will be read.
         */
        std::shared_ptr<Dataset> construct(const io::ExistingFile&, unsigned int expected_cols) override;
    };
}