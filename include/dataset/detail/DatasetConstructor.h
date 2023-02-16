#pragma once

#include <memory>
#include <string>

class Dataset;

namespace detail {
    /**
     * @brief Abstract base class for constructing datasets from files.
     */
    struct DatasetConstructor {
        /**
         * @brief Construct a dataset from a file.
         * 
         * @param path The path to the file.
         * @param expected_cols The expected number of columns. Any additional columns will be ignored.
         */
        virtual std::shared_ptr<Dataset> construct(std::string path, unsigned int expected_cols) = 0;
    };
}