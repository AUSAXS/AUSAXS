#pragma once

#include <dataset/detail/DatasetReader.h>

#include <unordered_set>

namespace detail {
    /**
     * @brief Contructs a dataset from a GROMACS XVG file.
     */
    struct XVGReader : DatasetReader {
        ~XVGReader() override = default;

        /**
         * @brief Construct a dataset from a GROMACS XVG file.
         * 
         * @param path The path to the file.
         * @param expected_cols The expected number of columns. Any additional columns will be ignored. If 0, all columns will be read.
         */
        std::unique_ptr<Dataset> construct(const io::ExistingFile&, unsigned int expected_cols) override;

        inline static std::unordered_set<std::string> extensions = {".xvg"};
    };
}