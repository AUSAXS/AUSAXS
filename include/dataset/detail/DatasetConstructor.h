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
         */
        virtual std::shared_ptr<Dataset> construct(std::string path) = 0;
    };
}