#pragma once

#include <dataset/detail/DatasetReader.h>

#include <unordered_set>

namespace ausaxs::detail {
    /**
     * @brief Contructs a dataset from a DAT file.
     *        Supports the formats x | y | yerr and x | y | yerr | xerr.
     */
    struct DATReader : DatasetReader {
        ~DATReader() override = default;

        /**
         * @brief Construct a dataset from a DAT file.
         * 
         * @param path The path to the file.
         * @param expected_cols The expected number of columns. Any additional columns will be ignored. If 0, all columns will be read.
         */
        std::unique_ptr<Dataset> construct(const io::ExistingFile&, unsigned int expected_cols) override;

        inline static std::unordered_set<std::string> extensions = {".dat", ".txt", ".rsr"};
    };
}