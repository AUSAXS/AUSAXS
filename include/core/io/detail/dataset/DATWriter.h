// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <dataset/DatasetFwd.h>
#include <io/IOFwd.h>

#include <string>

namespace ausaxs::io::detail::dataset {
    /**
     * @brief Writes a dataset to a DAT file.
     */
    struct DATWriter {
        /**
         * @brief Write a dataset to a DAT file.
         *
         * @param data The dataset to write.
         * @param path The path to the save location.
         * @param header The header for the file.
         */
        static void write(const Dataset& data, const io::File& path, const std::string& header);
    };
}
