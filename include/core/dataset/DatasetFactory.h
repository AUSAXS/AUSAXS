// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <dataset/Dataset.h>
#include <io/ExistingFile.h>

#include <memory>

namespace ausaxs::factory {
    /**
     * @brief Factory class for creating datasets from files.
     *        The factory will automatically determine the correct constructor to use based on the file extension.
     */
    struct DatasetFactory {
        /**
         * @brief Construct a dataset from a file.
         *        The returned dataset is guaranteed to be of the form: 
         *            [x | y] if @a filename is two-dimensional,
         *            [x | y | yerr] if @a filename is three-dimensional, or
         *            [x | y | yerr | xerr] if @a filename is four-dimensional.
         *        The return type depends on the @a expected_cols parameter. 
         *        If the file contains more columns than @a expected_cols, the extra columns will be ignored, and a warning will be printed.
         * 
         * @param file The path to the file.
         * @param expected_cols The expected number of columns. Any additional columns will be ignored. If 0, all columns will be read.
         */
        static std::unique_ptr<Dataset> construct(const io::ExistingFile& file, unsigned int expected_cols = 0);
    };
}