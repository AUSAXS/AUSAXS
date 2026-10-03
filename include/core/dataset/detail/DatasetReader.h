// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <dataset/DatasetFwd.h>
#include <io/ExistingFile.h>

#include <memory>
#include <string>

namespace ausaxs::detail {
    /**
     * @brief Abstract base class for constructing datasets from files.
     */
    struct DatasetReader {
        virtual ~DatasetReader() = default;

        /**
         * @brief Construct a dataset from a file.
         * 
         * @param path The path to the file.
         * @param expected_cols The expected number of columns. Any additional columns will be ignored.
         */
        virtual std::unique_ptr<Dataset> construct(const io::ExistingFile& path, int expected_cols) = 0;
    };

    /**
     * @brief The message for a data file with fewer columns than the caller needs.
     */
    inline std::string message_too_few_columns(const io::ExistingFile& path, int found, int expected) {
        std::string message = "\"" + path.str() + "\" has " + std::to_string(found) + " columns, but " + std::to_string(expected) + " are required";
        if (expected == 3) {
            message += " [q | I | Ierr]";
        }
        return message + ".";
    }
}