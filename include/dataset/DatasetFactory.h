#pragma once

#include <memory>

#include <dataset/Dataset.h>

namespace factory {
    /**
     * @brief Factory class for creating datasets from files.
     *        The factory will automatically determine the correct constructor to use based on the file extension.
     */
    struct DatasetFactory {
        /**
         * @brief Construct a dataset from a file.
         */
        static std::shared_ptr<Dataset> construct(std::string filename);
    };
}